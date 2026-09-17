







"""Parser for XML results returned by NCBI's Entrez Utilities.

This parser is used by the read() function in Bio.Entrez, and is not
intended be used directly.

The question is how to represent an XML file as Python objects. Some
XML files returned by NCBI look like lists, others look like dictionaries,
and others look like a mix of lists and dictionaries.

My approach is to classify each possible element in the XML as a plain
string, an integer, a list, a dictionary, or a structure. The latter is a
dictionary where the same key can occur multiple times; in Python, it is
represented as a dictionary where that key occurs once, pointing to a list
of values found in the XML file.

The parser then goes through the XML and creates the appropriate Python
object for each element. The different levels encountered in the XML are
preserved on the Python side. So a subelement of a subelement of an element
is a value in a dictionary that is stored in a list which is a value in
some other dictionary (or a value in a list which itself belongs to a list
which is a value in a dictionary, and so on). Attributes encountered in
the XML are stored as a dictionary in a member .attributes of each element,
and the tag name is saved in a member .tag.

To decide which kind of Python object corresponds to each element in the
XML, the parser analyzes the DTD referred at the top of (almost) every
XML file returned by the Entrez Utilities. This is preferred over a hand-
written solution, since the number of DTDs is rather large and their
contents may change over time. About half the code in this parser deals
with parsing the DTD, and the other half with the XML itself.
"""

import os
import warnings
import xml.etree.ElementTree as ET
from collections import Counter
from io import BytesIO
from urllib.parse import urlparse
from urllib.request import urlopen
from xml.parsers import expat
from xml.sax.saxutils import escape

from Bio import StreamModeError





class NoneElement:
    """NCBI Entrez XML element mapped to None."""

    def __init__(self, tag, attributes, key):
        """Create a NoneElement."""
        self.tag = tag
        self.key = key
        self.attributes = attributes

    def __repr__(self):
        """Return a string representation of the object."""
        try:
            attributes = self.attributes
        except AttributeError:
            return "NoneElement"
        return "NoneElement(attributes=%r)" % attributes

    def __eq__(self, other):
        if isinstance(other, NoneElement):
            return True
        return False


class IntegerElement(int):
    """NCBI Entrez XML element mapped to an integer."""

    def __new__(cls, value, *args, **kwargs):
        """Create an IntegerElement."""
        return int.__new__(cls, value)

    def __init__(self, value, tag, attributes, key):
        """Initialize an IntegerElement."""
        self.tag = tag
        self.attributes = attributes
        self.key = key

    def __repr__(self):
        """Return a string representation of the object."""
        text = int.__repr__(self)
        try:
            attributes = self.attributes
        except AttributeError:
            return text
        return f"IntegerElement({text}, attributes={attributes!r})"


class StringElement(str):
    """NCBI Entrez XML element mapped to a string."""

    def __new__(cls, value, *args, **kwargs):
        """Create a StringElement."""
        return str.__new__(cls, value)

    def __init__(self, value, tag, attributes, key):
        """Initialize a StringElement."""
        self.tag = tag
        self.attributes = attributes
        self.key = key

    def __repr__(self):
        """Return a string representation of the object."""
        text = str.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return f"StringElement({text}, attributes={attributes!r})"


class ListElement(list):
    """NCBI Entrez XML element mapped to a list."""

    def __init__(self, tag, attributes, allowed_tags, key=None):
        """Create a ListElement."""
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes
        self.allowed_tags = allowed_tags

    def __repr__(self):
        """Return a string representation of the object."""
        text = list.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return f"ListElement({text}, attributes={attributes!r})"

    def store(self, value):
        """Append an element to the list, checking tags."""
        key = value.key
        if self.allowed_tags is not None and key not in self.allowed_tags:
            raise ValueError("Unexpected item '%s' in list" % key)
        del value.key
        self.append(value)


class DictionaryElement(dict):
    """NCBI Entrez XML element mapped to a dictionaray."""

    def __init__(self, tag, attrs, allowed_tags, repeated_tags=None, key=None):
        """Create a DictionaryElement."""
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attrs
        self.allowed_tags = allowed_tags
        self.repeated_tags = repeated_tags
        if repeated_tags:
            for key in repeated_tags:
                self[key] = []

    def __repr__(self):
        """Return a string representation of the object."""
        text = dict.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return f"DictElement({text}, attributes={attributes!r})"

    def store(self, value):
        """Add an entry to the dictionary, checking tags."""
        key = value.key
        tag = value.tag
        if self.allowed_tags is not None and tag not in self.allowed_tags:
            raise ValueError("Unexpected item '%s' in dictionary" % key)
        del value.key
        if self.repeated_tags and key in self.repeated_tags:
            self[key].append(value)
        else:
            self[key] = value


class OrderedListElement(list):
    """NCBI Entrez XML element mapped to a list of lists.

    OrderedListElement is used to describe a list of repeating elements such as
    A, B, C, A, B, C, A, B, C ... where each set of A, B, C forms a group. This
    is then stored as [[A, B, C], [A, B, C], [A, B, C], ...]
    """

    def __init__(self, tag, attributes, allowed_tags, first_tag, key=None):
        """Create an OrderedListElement."""
        self.tag = tag
        if key is None:
            self.key = tag
        else:
            self.key = key
        self.attributes = attributes
        self.allowed_tags = allowed_tags
        self.first_tag = first_tag

    def __repr__(self):
        """Return a string representation of the object."""
        text = list.__repr__(self)
        attributes = self.attributes
        if not attributes:
            return text
        return f"OrderedListElement({text}, attributes={attributes!r})"

    def store(self, value):
        """Append an element to the list, checking tags."""
        key = value.key
        if self.allowed_tags is not None and key not in self.allowed_tags:
            raise ValueError("Unexpected item '%s' in list" % key)
        if key == self.first_tag:
            self.append([])
        self[-1].append(value)


class ErrorElement(str):
    """NCBI Entrez XML element containing an error message."""

    def __new__(cls, value, *args, **kwargs):
        """Create an ErrorElement."""
        return str.__new__(cls, value)

    def __init__(self, value, tag):
        """Initialize an ErrorElement."""
        self.tag = tag
        self.key = tag

    def __repr__(self):
        """Return the error message as a string."""
        text = str.__repr__(self)
        return f"ErrorElement({text})"


class NotXMLError(ValueError):
    """Failed to parse file as XML."""

    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to parse the XML data (%s). Please make sure that the input data "
            "are in XML format." % self.msg
        )


class CorruptedXMLError(ValueError):
    """Corrupted XML."""

    def __init__(self, message):
        """Initialize the class."""
        self.msg = message

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to parse the XML data (%s). Please make sure that the input data "
            "are not corrupted." % self.msg
        )


class ValidationError(ValueError):
    """XML tag found which was not defined in the DTD.

    Validating parsers raise this error if the parser finds a tag in the XML
    that is not defined in the DTD. Non-validating parsers do not raise this
    error. The Bio.Entrez.read and Bio.Entrez.parse functions use validating
    parsers by default (see those functions for more information).
    """

    def __init__(self, name):
        """Initialize the class."""
        self.name = name

    def __str__(self):
        """Return a string summary of the exception."""
        return (
            "Failed to find tag '%s' in the DTD. To skip all tags that "
            "are not represented in the DTD, please call Bio.Entrez.read "
            "or Bio.Entrez.parse with validate=False." % self.name
        )


class DataHandlerMeta(type):
    """A metaclass is needed until Python supports @classproperty."""

    def __init__(cls, *args, **kwargs):
        """Initialize the class."""
        from Bio import Entrez

        try:
            cls.directory = Entrez.local_cache  
        except PermissionError:
            cls._directory = Entrez.local_cache  
        del Entrez

    @property
    def directory(cls):
        """Directory for caching XSD and DTD files."""
        return cls._directory

    @directory.setter
    def directory(cls, value):
        """Set a custom directory for the local DTD/XSD directories."""
        if value is None:
            import platform

            if platform.system() == "Windows":
                value = os.path.join(os.getenv("APPDATA"), "biopython")
            else:  
                home = os.path.expanduser("~")
                value = os.path.join(home, ".config", "biopython")
        
        cls.local_dtd_dir = os.path.join(value, "Bio", "Entrez", "DTDs")
        os.makedirs(cls.local_dtd_dir, exist_ok=True)
        
        cls.local_xsd_dir = os.path.join(value, "Bio", "Entrez", "XSDs")
        os.makedirs(cls.local_xsd_dir, exist_ok=True)
        
        
        cls._directory = value


class DataHandler(metaclass=DataHandlerMeta):
    """Data handler for parsing NCBI XML from Entrez."""

    from Bio import Entrez

    global_dtd_dir = os.path.join(Entrez.__path__[0], "DTDs")
    global_xsd_dir = os.path.join(Entrez.__path__[0], "XSDs")
    local_dtd_dir = None
    local_xsd_dir = None

    del Entrez

    def __init__(self, validate, escape, ignore_errors):
        """Create a DataHandler object."""
        self.dtd_urls = []
        self.element = None
        self.level = 0
        self.bypass_url_security = False
        self.data = []
        self.attributes = None
        self.allowed_tags = None
        self.constructors = {}
        self.strings = {}
        self.items = set()
        self.errors = set()
        self.validating = validate
        self.ignore_errors = ignore_errors
        self.parser = expat.ParserCreate(namespace_separator=" ")
        self.parser.SetParamEntityParsing(expat.XML_PARAM_ENTITY_PARSING_ALWAYS)
        self.parser.XmlDeclHandler = self.xmlDeclHandler
        self.schema_namespace = None
        self.namespace_level = Counter()
        self.namespace_prefix = {}
        if escape:
            self.characterDataHandler = self.characterDataHandlerEscape
        else:
            self.characterDataHandler = self.characterDataHandlerRaw

    def read(self, source):
        """Set up the parser and let it read the XML results."""
        
        
        try:
            stream = open(source, "rb")
        except TypeError:  
            if source.read(0) != b"":
                raise StreamModeError(
                    "the XML file must be opened in binary mode."
                ) from None
            stream = source
        if stream.read(0) != b"":
            raise TypeError("file should be opened in binary mode")
        try:
            self.parser.ParseFile(stream)
        except expat.ExpatError as e:
            if self.parser.StartElementHandler:
                
                
                
                raise CorruptedXMLError(e) from None
            else:
                
                
                raise NotXMLError(e) from None
        finally:
            if stream is not source:
                stream.close()
        try:
            record = self.record
        except AttributeError:
            if self.parser.StartElementHandler:
                
                
                
                raise RuntimeError(
                    "Failed to parse the XML file correctly, possibly due to a bug "
                    "in Bio.Entrez. Please contact the Biopython developers via "
                    "the mailing list or GitHub for assistance."
                ) from None
            else:
                
                
                raise NotXMLError("XML declaration not found") from None
        else:
            del record.key
            return record

    def parse(self, source):
        """Set up the parser and let it read the XML results."""
        
        
        
        
        
        
        
        
        
        
        
        
        
        try:
            stream = open(source, "rb")
        except TypeError:  
            if source.read(0) != b"":
                raise StreamModeError(
                    "the XML file must be opened in binary mode."
                ) from None
            stream = source
        if stream.read(0) != b"":
            raise TypeError("file should be opened in binary mode")
        BLOCK = 1024
        try:
            while True:
                
                data = stream.read(BLOCK)
                self.parser.Parse(data, False)
                try:
                    records = self.record
                except AttributeError:
                    if self.parser.StartElementHandler:
                        
                        
                        

                        raise RuntimeError(
                            "Failed to parse the XML file correctly, possibly due to a "
                            "bug in Bio.Entrez. Please contact the Biopython "
                            "developers via the mailing list or GitHub for assistance."
                        ) from None
                    else:
                        
                        
                        raise NotXMLError("XML declaration not found") from None

                if not isinstance(records, list):
                    raise ValueError(
                        "The XML file does not represent a list. Please use "
                        "Entrez.read instead of Entrez.parse."
                    )

                if not data:
                    break

                while len(records) >= 2:
                    
                    
                    record = records.pop(0)
                    yield record

        except expat.ExpatError as e:
            if self.parser.StartElementHandler:
                
                
                
                raise CorruptedXMLError(e) from None
            else:
                
                
                raise NotXMLError(e) from None
        finally:
            if stream is not source:
                stream.close()

        
        self.parser = None
        if self.element is not None:
            
            raise CorruptedXMLError("Premature end of data")

        
        yield from records

    def xmlDeclHandler(self, version, encoding, standalone):
        """Set XML handlers when an XML declaration is found."""
        self.parser.CharacterDataHandler = self.characterDataHandler
        self.parser.ExternalEntityRefHandler = self.externalEntityRefHandler
        self.parser.StartNamespaceDeclHandler = self.startNamespaceDeclHandler
        self.parser.EndNamespaceDeclHandler = self.endNamespaceDeclHandler
        self.parser.StartElementHandler = self.handleMissingDocumentDefinition

    def handleMissingDocumentDefinition(self, tag, attrs):
        """Raise an Exception if neither a DTD nor an XML Schema is found."""
        raise ValueError(
            "As the XML data contained neither a Document Type Definition (DTD) nor an XML Schema, Bio.Entrez is unable to parse these data. We recommend using a generic XML parser from the Python standard library instead, for example ElementTree."
        )

    def startNamespaceDeclHandler(self, prefix, uri):
        """Handle start of an XML namespace declaration."""
        if prefix == "xsi":
            
            self.schema_namespace = uri
            self.parser.StartElementHandler = self.schemaHandler
        else:
            
            
            
            
            
            
            
            if prefix == "mml":
                assert uri == "http://www.w3.org/1998/Math/MathML"
            elif prefix == "xlink":
                assert uri == "http://www.w3.org/1999/xlink"
            elif prefix == "ali":
                assert uri.rstrip("/") == "http://www.niso.org/schemas/ali/1.0"
            else:
                raise ValueError(f"Unknown prefix '{prefix}' with uri '{uri}'")
            self.namespace_level[prefix] += 1
            self.namespace_prefix[uri] = prefix

    def endNamespaceDeclHandler(self, prefix):
        """Handle end of an XML namespace declaration."""
        if prefix != "xsi":
            self.namespace_level[prefix] -= 1
            if self.namespace_level[prefix] == 0:
                for key, value in self.namespace_prefix.items():
                    if value == prefix:
                        break
                else:
                    raise RuntimeError("Failed to find namespace prefix")
                del self.namespace_prefix[key]

    def schemaHandler(self, name, attrs):
        """Process the XML schema (before processing the element)."""
        key = "%s noNamespaceSchemaLocation" % self.schema_namespace
        schema = attrs[key]
        self.verify_security(schema)
        handle = self.open_xsd_file(os.path.basename(schema))
        
        if not handle:
            handle = urlopen(schema)
            text = handle.read()
            self.save_xsd_file(os.path.basename(schema), text)
            handle.close()
            self.parse_xsd(ET.fromstring(text))
        else:
            self.parse_xsd(ET.fromstring(handle.read()))
            handle.close()
        
        self.startElementHandler(name, attrs)
        
        self.parser.StartElementHandler = self.startElementHandler

    def startElementHandler(self, tag, attrs):
        """Handle start of an XML element."""
        prefix = None
        if self.namespace_prefix:
            try:
                uri, name = tag.split()
            except ValueError:
                pass
            else:
                prefix = self.namespace_prefix[uri]
                tag = f"{prefix}:{name}"
        if tag in self.items:
            assert tag == "Item"
            name = attrs["Name"]
            itemtype = attrs["Type"]
            del attrs["Type"]
            if itemtype == "Structure":
                del attrs["Name"]
                element = DictionaryElement(
                    name, attrs, allowed_tags=None, repeated_tags=None
                )
                parent = self.element
                element.parent = parent
                
                if parent is None:
                    self.record = element
                else:
                    parent.store(element)
                self.element = element
                self.parser.EndElementHandler = self.endElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            elif name in ("ArticleIds", "History"):
                del attrs["Name"]
                allowed_tags = None  
                repeated_tags = frozenset(["pubmed", "medline"])
                element = DictionaryElement(
                    tag,
                    attrs,
                    allowed_tags=allowed_tags,
                    repeated_tags=repeated_tags,
                    key=name,
                )
                parent = self.element
                element.parent = parent
                
                if parent is None:
                    self.record = element
                else:
                    parent.store(element)
                self.element = element
                self.parser.EndElementHandler = self.endElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            elif itemtype == "List":
                del attrs["Name"]
                allowed_tags = None  
                element = ListElement(tag, attrs, allowed_tags, name)
                parent = self.element
                element.parent = parent
                if self.element is None:
                    
                    self.record = element
                else:
                    parent.store(element)
                self.element = element
                self.parser.EndElementHandler = self.endElementHandler
                self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            elif itemtype == "Integer":
                self.parser.EndElementHandler = self.endIntegerElementHandler
                self.parser.CharacterDataHandler = self.characterDataHandler
                self.attributes = attrs
            elif itemtype in ("String", "Unknown", "Date", "Enumerator"):
                assert self.attributes is None
                self.attributes = attrs
                self.parser.StartElementHandler = self.startRawElementHandler
                self.parser.EndElementHandler = self.endStringElementHandler
                self.parser.CharacterDataHandler = self.characterDataHandler
            else:
                raise ValueError("Unknown item type %s" % name)
        elif tag in self.errors:
            self.parser.EndElementHandler = self.endErrorElementHandler
            self.parser.CharacterDataHandler = self.characterDataHandler
        elif tag in self.strings:
            self.parser.StartElementHandler = self.startRawElementHandler
            self.parser.EndElementHandler = self.endStringElementHandler
            self.parser.CharacterDataHandler = self.characterDataHandler
            assert self.allowed_tags is None
            self.allowed_tags = self.strings[tag]
            assert self.attributes is None
            self.attributes = attrs
        elif tag in self.constructors:
            cls, allowed_tags = self.constructors[tag]
            element = cls(tag, attrs, *allowed_tags)
            parent = self.element
            element.parent = parent
            if parent is None:
                
                self.record = element
            else:
                parent.store(element)
            self.element = element
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
        else:
            
            if tag == "processing-meta":
                terms = []
                dtd_version = "1.3"
                if attrs["tagset-family"] == "jats":
                    terms.append("JATS")
                if attrs["base-tagset"] == "archiving":
                    term = "archivearticle" + dtd_version.replace(".", "-")
                    terms.append(term)
                if attrs.get("mathml-version") == "3.0":
                    terms.append("mathml3")
                basename = "-".join(terms)
                url = f"https://{attrs['tagset-family']}.nlm.nih.gov/{attrs['base-tagset']}/{dtd_version}/{basename}.dtd"
                self.xmlDeclHandler(None, None, None)
                self.externalEntityRefHandler(None, None, url, None)
                
            elif self.validating:
                raise ValidationError(tag)
            
            self.parser.StartElementHandler = self.startSkipElementHandler
            self.parser.EndElementHandler = self.endSkipElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            self.level = 1

    def startRawElementHandler(self, name, attrs):
        """Handle start of an XML raw element."""
        
        prefix = None
        if self.namespace_prefix:
            try:
                uri, name = name.split()
            except ValueError:
                pass
            else:
                prefix = self.namespace_prefix[uri]
                if self.namespace_level[prefix] == 1:
                    attrs = {"xmlns": uri}
        if prefix:
            key = f"{prefix}:{name}"
        else:
            key = name
        
        
        tag = "<%s" % name
        for key, value in attrs.items():
            tag += f' {key}="{value}"'
        tag += ">"
        self.data.append(tag)
        self.parser.EndElementHandler = self.endRawElementHandler
        self.level += 1

    def startSkipElementHandler(self, name, attrs):
        """Handle start of an XML skip element."""
        self.level += 1

    def endStringElementHandler(self, tag):
        """Handle end of an XML string element."""
        element = self.element
        if element is not None:
            self.parser.StartElementHandler = self.startElementHandler
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
        data = "".join(self.data)
        self.data = []
        attributes = self.attributes
        self.attributes = None
        if self.namespace_prefix:
            try:
                uri, name = tag.split()
            except ValueError:
                pass
            else:
                prefix = self.namespace_prefix[uri]
                tag = f"{prefix}:{name}"
        if tag in self.items:
            assert tag == "Item"
            key = attributes["Name"]
            del attributes["Name"]
        else:
            key = tag
        value = StringElement(data, tag, attributes, key)
        if element is None:
            self.record = element
        else:
            element.store(value)
        self.allowed_tags = None

    def endRawElementHandler(self, name):
        """Handle end of an XML raw element."""
        self.level -= 1
        if self.level == 0:
            self.parser.EndElementHandler = self.endStringElementHandler
        if self.namespace_prefix:
            try:
                uri, name = name.split()
            except ValueError:
                pass
        tag = "</%s>" % name
        self.data.append(tag)

    def endSkipElementHandler(self, name):
        """Handle end of an XML skip element."""
        self.level -= 1
        if self.level == 0:
            self.parser.StartElementHandler = self.startElementHandler
            self.parser.EndElementHandler = self.endElementHandler

    def endErrorElementHandler(self, tag):
        """Handle end of an XML error element."""
        element = self.element
        if element is not None:
            self.parser.StartElementHandler = self.startElementHandler
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
        data = "".join(self.data)
        if data == "":
            return
        if self.ignore_errors is False:
            raise RuntimeError(data)
        self.data = []
        value = ErrorElement(data, tag)
        if element is None:
            self.record = element
        else:
            element.store(value)

    def endElementHandler(self, name):
        """Handle end of an XML element."""
        element = self.element
        self.element = element.parent
        del element.parent

    def endIntegerElementHandler(self, tag):
        """Handle end of an XML integer element."""
        attributes = self.attributes
        self.attributes = None
        assert tag == "Item"
        key = attributes["Name"]
        del attributes["Name"]
        if self.data:
            value = int("".join(self.data))
            self.data = []
            value = IntegerElement(value, tag, attributes, key)
        else:
            value = NoneElement(tag, attributes, key)
        element = self.element
        if element is None:
            self.record = value
        else:
            self.parser.EndElementHandler = self.endElementHandler
            self.parser.CharacterDataHandler = self.skipCharacterDataHandler
            if value is None:
                return
            element.store(value)

    def characterDataHandlerRaw(self, content):
        """Handle character data as-is (raw)."""
        self.data.append(content)

    def characterDataHandlerEscape(self, content):
        """Handle character data by encoding it."""
        content = escape(content)
        self.data.append(content)

    def skipCharacterDataHandler(self, content):
        """Handle character data by skipping it."""

    def parse_xsd(self, root):
        """Parse an XSD file."""
        prefix = "{http://www.w3.org/2001/XMLSchema}"
        for element in root:
            isSimpleContent = False
            attribute_keys = []
            keys = []
            multiple = []
            assert element.tag == prefix + "element"
            name = element.attrib["name"]
            assert len(element) == 1
            complexType = element[0]
            assert complexType.tag == prefix + "complexType"
            for component in complexType:
                tag = component.tag
                if tag == prefix + "attribute":
                    
                    attribute_keys.append(component.attrib["name"])
                elif tag == prefix + "sequence":
                    maxOccurs = component.attrib.get("maxOccurs", "1")
                    for key in component:
                        assert key.tag == prefix + "element"
                        ref = key.attrib["ref"]
                        keys.append(ref)
                        if maxOccurs != "1" or key.attrib.get("maxOccurs", "1") != "1":
                            multiple.append(ref)
                elif tag == prefix + "simpleContent":
                    assert len(component) == 1
                    extension = component[0]
                    assert extension.tag == prefix + "extension"
                    assert extension.attrib["base"] == "xs:string"
                    for attribute in extension:
                        assert attribute.tag == prefix + "attribute"
                        
                        attribute_keys.append(attribute.attrib["name"])
                    isSimpleContent = True
            allowed_tags = frozenset(keys)
            if len(keys) == 1 and keys == multiple:
                assert not isSimpleContent
                args = (allowed_tags,)
                self.constructors[name] = (ListElement, args)
            elif len(keys) >= 1:
                assert not isSimpleContent
                repeated_tags = frozenset(multiple)
                args = (allowed_tags, repeated_tags)
                self.constructors[name] = (DictionaryElement, args)
            else:
                self.strings[name] = allowed_tags

    def elementDecl(self, name, model):
        """Call a call-back function for each element declaration in a DTD.

        This is used for each element declaration in a DTD like::

            <!ELEMENT       name          (...)>

        The purpose of this function is to determine whether this element
        should be regarded as a string, integer, list, dictionary, structure,
        or error.
        """
        if name.upper() == "ERROR":
            self.errors.add(name)
            return
        if name == "Item" and model == (
            expat.model.XML_CTYPE_MIXED,
            expat.model.XML_CQUANT_REP,
            None,
            ((expat.model.XML_CTYPE_NAME, expat.model.XML_CQUANT_NONE, "Item", ()),),
        ):
            
            
            self.items.add(name)
            return
        
        while (
            model[0] in (expat.model.XML_CTYPE_SEQ, expat.model.XML_CTYPE_CHOICE)
            and model[1] in (expat.model.XML_CQUANT_NONE, expat.model.XML_CQUANT_OPT)
            and len(model[3]) == 1
        ):
            model = model[3][0]
        
        if model[0] in (expat.model.XML_CTYPE_MIXED, expat.model.XML_CTYPE_EMPTY):
            if model[1] == expat.model.XML_CQUANT_REP:
                children = model[3]
                allowed_tags = frozenset(child[2] for child in children)
            else:
                allowed_tags = frozenset()
            self.strings[name] = allowed_tags
            return
        
        if model == (expat.model.XML_CTYPE_ANY, expat.model.XML_CQUANT_NONE, None, ()):
            allowed_tags = None
            repeated_tags = None
            args = (allowed_tags, repeated_tags)
            self.constructors[name] = (DictionaryElement, args)
            return
        
        if model[0] in (
            expat.model.XML_CTYPE_CHOICE,
            expat.model.XML_CTYPE_SEQ,
        ) and model[1] in (expat.model.XML_CQUANT_PLUS, expat.model.XML_CQUANT_REP):
            children = model[3]
            allowed_tags = frozenset(child[2] for child in children)
            if model[0] == expat.model.XML_CTYPE_SEQ:
                if len(children) > 1:
                    assert model[1] == expat.model.XML_CQUANT_PLUS
                    first_child = children[0]
                    assert first_child[1] == expat.model.XML_CQUANT_NONE
                    first_tag = first_child[2]
                    args = allowed_tags, first_tag
                    self.constructors[name] = (OrderedListElement, args)
                    return
                assert len(children) == 1
            self.constructors[name] = (ListElement, (allowed_tags,))
            return
        
        
        
        
        
        
        
        single = []
        multiple = []
        errors = []
        
        

        def count(model):
            quantifier, key, children = model[1:]
            if key is None:
                if quantifier in (
                    expat.model.XML_CQUANT_PLUS,
                    expat.model.XML_CQUANT_REP,
                ):
                    for child in children:
                        multiple.append(child[2])
                else:
                    for child in children:
                        count(child)
            elif key.upper() == "ERROR":
                errors.append(key)
            else:
                if quantifier in (
                    expat.model.XML_CQUANT_NONE,
                    expat.model.XML_CQUANT_OPT,
                ):
                    single.append(key)
                elif quantifier in (
                    expat.model.XML_CQUANT_PLUS,
                    expat.model.XML_CQUANT_REP,
                ):
                    multiple.append(key)

        count(model)
        if len(single) == 0 and len(multiple) == 1:
            allowed_tags = frozenset(multiple + errors)
            self.constructors[name] = (ListElement, (allowed_tags,))
        else:
            allowed_tags = frozenset(single + multiple + errors)
            repeated_tags = frozenset(multiple)
            args = (allowed_tags, repeated_tags)
            self.constructors[name] = (DictionaryElement, args)

    def open_dtd_file(self, filename):
        """Open specified DTD file."""
        if DataHandler.local_dtd_dir is not None:
            path = os.path.join(DataHandler.local_dtd_dir, filename)
            try:
                handle = open(path, "rb")
            except FileNotFoundError:
                pass
            else:
                return handle
        path = os.path.join(DataHandler.global_dtd_dir, filename)
        try:
            handle = open(path, "rb")
        except FileNotFoundError:
            pass
        else:
            return handle
        return None

    def open_xsd_file(self, filename):
        """Open specified XSD file."""
        if DataHandler.local_xsd_dir is not None:
            path = os.path.join(DataHandler.local_xsd_dir, filename)
            try:
                handle = open(path, "rb")
            except FileNotFoundError:
                pass
            else:
                return handle
        path = os.path.join(DataHandler.global_xsd_dir, filename)
        try:
            handle = open(path, "rb")
        except FileNotFoundError:
            pass
        else:
            return handle
        return None

    def save_dtd_file(self, filename, text):
        """Save DTD file to cache."""
        if DataHandler.local_dtd_dir is None:
            return
        path = os.path.join(DataHandler.local_dtd_dir, filename)
        try:
            handle = open(path, "wb")
        except OSError:
            warnings.warn(f"Failed to save {filename} at {path}")
        else:
            handle.write(text)
            handle.close()

    def save_xsd_file(self, filename, text):
        """Save XSD file to cache."""
        if DataHandler.local_xsd_dir is None:
            return
        path = os.path.join(DataHandler.local_xsd_dir, filename)
        try:
            handle = open(path, "wb")
        except OSError:
            warnings.warn(f"Failed to save {filename} at {path}")
        else:
            handle.write(text)
            handle.close()

    def verify_security(self, url, verify_hostname=True):
        """Check if the given URL is from a trustable source.

        When ``self.bypass_url_security`` evaluates to ``True``,
        all URL security checks will be skipped.
        """
        if not self.bypass_url_security:
            parts = urlparse(url)
            scheme = parts.scheme
            hostname = parts.hostname
            if scheme != "https" or (
                verify_hostname and not hostname.endswith(".nlm.nih.gov")
            ):
                raise ValueError(f"Expected secure URL to NCBI, found {url!r}")

    def externalEntityRefHandler(self, context, base, systemId, publicId):
        """Handle external entity reference in order to cache DTD locally.

        The purpose of this function is to load the DTD locally, instead
        of downloading it from the URL specified in the XML. Using the local
        DTD results in much faster parsing. If the DTD is not found locally,
        we try to download it. If new DTDs become available from NCBI,
        putting them in Bio/Entrez/DTDs will allow the parser to see them.
        """
        urlinfo = urlparse(systemId)
        if urlinfo.scheme in ["http", "https", "ftp"]:
            
            url = systemId
        elif urlinfo.scheme == "":
            
            
            try:
                source = self.dtd_urls[-1]
            except IndexError:
                
                
                source = "https://www.ncbi.nlm.nih.gov/dtd/"
            else:
                source = os.path.dirname(source)
            
            url = source.rstrip("/") + "/" + systemId
        else:
            raise ValueError("Unexpected URL scheme %r" % urlinfo.scheme)

        
        
        
        
        
        self.verify_security(url, verify_hostname=not self.dtd_urls)

        
        
        
        self.dtd_urls.append(url)
        try:
            
            location, filename = os.path.split(systemId)
            handle = self.open_dtd_file(filename)
            if not handle:
                
                
                try:
                    handle = urlopen(url)
                except OSError:
                    raise RuntimeError(
                        f"Failed to access {filename} at {url}"
                    ) from None
                text = handle.read()
                handle.close()
                self.save_dtd_file(filename, text)
                handle = BytesIO(text)

            parser = self.parser.ExternalEntityParserCreate(context)
            parser.ElementDeclHandler = self.elementDecl
            parser.ParseFile(handle)
            handle.close()
        finally:
            self.dtd_urls.pop()

        self.parser.StartElementHandler = self.startElementHandler
        return 1
