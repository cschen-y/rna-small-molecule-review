#define PY_SSIZE_T_CLEAN
#include "Python.h"


static const char bases[][4] = {"TTTT",                   
                                "TTTC",                   
                                "TTTA",                   
                                "TTTG",                   
                                "TTCT",                   
                                "TTCC",                   
                                "TTCA",                   
                                "TTCG",                   
                                "TTAT",                   
                                "TTAC",                   
                                "TTAA",                   
                                "TTAG",                   
                                "TTGT",                   
                                "TTGC",                   
                                "TTGA",                   
                                "TTGG",                   
                                "TCTT",                   
                                "TCTC",                   
                                "TCTA",                   
                                "TCTG",                   
                                "TCCT",                   
                                "TCCC",                   
                                "TCCA",                   
                                "TCCG",                   
                                "TCAT",                   
                                "TCAC",                   
                                "TCAA",                   
                                "TCAG",                   
                                "TCGT",                   
                                "TCGC",                   
                                "TCGA",                   
                                "TCGG",                   
                                "TATT",                   
                                "TATC",                   
                                "TATA",                   
                                "TATG",                   
                                "TACT",                   
                                "TACC",                   
                                "TACA",                   
                                "TACG",                   
                                "TAAT",                   
                                "TAAC",                   
                                "TAAA",                   
                                "TAAG",                   
                                "TAGT",                   
                                "TAGC",                   
                                "TAGA",                   
                                "TAGG",                   
                                "TGTT",                   
                                "TGTC",                   
                                "TGTA",                   
                                "TGTG",                   
                                "TGCT",                   
                                "TGCC",                   
                                "TGCA",                   
                                "TGCG",                   
                                "TGAT",                   
                                "TGAC",                   
                                "TGAA",                   
                                "TGAG",                   
                                "TGGT",                   
                                "TGGC",                   
                                "TGGA",                   
                                "TGGG",                   
                                "CTTT",                   
                                "CTTC",                   
                                "CTTA",                   
                                "CTTG",                   
                                "CTCT",                   
                                "CTCC",                   
                                "CTCA",                   
                                "CTCG",                   
                                "CTAT",                   
                                "CTAC",                   
                                "CTAA",                   
                                "CTAG",                   
                                "CTGT",                   
                                "CTGC",                   
                                "CTGA",                   
                                "CTGG",                   
                                "CCTT",                   
                                "CCTC",                   
                                "CCTA",                   
                                "CCTG",                   
                                "CCCT",                   
                                "CCCC",                   
                                "CCCA",                   
                                "CCCG",                   
                                "CCAT",                   
                                "CCAC",                   
                                "CCAA",                   
                                "CCAG",                   
                                "CCGT",                   
                                "CCGC",                   
                                "CCGA",                   
                                "CCGG",                   
                                "CATT",                   
                                "CATC",                   
                                "CATA",                   
                                "CATG",                   
                                "CACT",                   
                                "CACC",                   
                                "CACA",                   
                                "CACG",                   
                                "CAAT",                   
                                "CAAC",                   
                                "CAAA",                   
                                "CAAG",                   
                                "CAGT",                   
                                "CAGC",                   
                                "CAGA",                   
                                "CAGG",                   
                                "CGTT",                   
                                "CGTC",                   
                                "CGTA",                   
                                "CGTG",                   
                                "CGCT",                   
                                "CGCC",                   
                                "CGCA",                   
                                "CGCG",                   
                                "CGAT",                   
                                "CGAC",                   
                                "CGAA",                   
                                "CGAG",                   
                                "CGGT",                   
                                "CGGC",                   
                                "CGGA",                   
                                "CGGG",                   
                                "ATTT",                   
                                "ATTC",                   
                                "ATTA",                   
                                "ATTG",                   
                                "ATCT",                   
                                "ATCC",                   
                                "ATCA",                   
                                "ATCG",                   
                                "ATAT",                   
                                "ATAC",                   
                                "ATAA",                   
                                "ATAG",                   
                                "ATGT",                   
                                "ATGC",                   
                                "ATGA",                   
                                "ATGG",                   
                                "ACTT",                   
                                "ACTC",                   
                                "ACTA",                   
                                "ACTG",                   
                                "ACCT",                   
                                "ACCC",                   
                                "ACCA",                   
                                "ACCG",                   
                                "ACAT",                   
                                "ACAC",                   
                                "ACAA",                   
                                "ACAG",                   
                                "ACGT",                   
                                "ACGC",                   
                                "ACGA",                   
                                "ACGG",                   
                                "AATT",                   
                                "AATC",                   
                                "AATA",                   
                                "AATG",                   
                                "AACT",                   
                                "AACC",                   
                                "AACA",                   
                                "AACG",                   
                                "AAAT",                   
                                "AAAC",                   
                                "AAAA",                   
                                "AAAG",                   
                                "AAGT",                   
                                "AAGC",                   
                                "AAGA",                   
                                "AAGG",                   
                                "AGTT",                   
                                "AGTC",                   
                                "AGTA",                   
                                "AGTG",                   
                                "AGCT",                   
                                "AGCC",                   
                                "AGCA",                   
                                "AGCG",                   
                                "AGAT",                   
                                "AGAC",                   
                                "AGAA",                   
                                "AGAG",                   
                                "AGGT",                   
                                "AGGC",                   
                                "AGGA",                   
                                "AGGG",                   
                                "GTTT",                   
                                "GTTC",                   
                                "GTTA",                   
                                "GTTG",                   
                                "GTCT",                   
                                "GTCC",                   
                                "GTCA",                   
                                "GTCG",                   
                                "GTAT",                   
                                "GTAC",                   
                                "GTAA",                   
                                "GTAG",                   
                                "GTGT",                   
                                "GTGC",                   
                                "GTGA",                   
                                "GTGG",                   
                                "GCTT",                   
                                "GCTC",                   
                                "GCTA",                   
                                "GCTG",                   
                                "GCCT",                   
                                "GCCC",                   
                                "GCCA",                   
                                "GCCG",                   
                                "GCAT",                   
                                "GCAC",                   
                                "GCAA",                   
                                "GCAG",                   
                                "GCGT",                   
                                "GCGC",                   
                                "GCGA",                   
                                "GCGG",                   
                                "GATT",                   
                                "GATC",                   
                                "GATA",                   
                                "GATG",                   
                                "GACT",                   
                                "GACC",                   
                                "GACA",                   
                                "GACG",                   
                                "GAAT",                   
                                "GAAC",                   
                                "GAAA",                   
                                "GAAG",                   
                                "GAGT",                   
                                "GAGC",                   
                                "GAGA",                   
                                "GAGG",                   
                                "GGTT",                   
                                "GGTC",                   
                                "GGTA",                   
                                "GGTG",                   
                                "GGCT",                   
                                "GGCC",                   
                                "GGCA",                   
                                "GGCG",                   
                                "GGAT",                   
                                "GGAC",                   
                                "GGAA",                   
                                "GGAG",                   
                                "GGGT",                   
                                "GGGC",                   
                                "GGGA",                   
                                "GGGG",                   
                               };

static int
extract(const unsigned char* bytes, Py_ssize_t byteSize,
        Py_ssize_t start, Py_ssize_t end, char sequence[]) {
    Py_ssize_t i;
    const Py_ssize_t size = end - start;
    const Py_ssize_t byteStart = start / 4;
    const Py_ssize_t byteEnd = (end + 3) / 4;

    if (byteSize != byteEnd - byteStart) {
        PyErr_Format(PyExc_RuntimeError,
                     "unexpected number of bytes %u (expected %u)",
                     byteSize, byteEnd - byteStart);
        return -1;
    }

    start -= byteStart * 4;
    if (byteStart + 1 == byteEnd) {
                           
        memcpy(sequence, &(bases[*bytes][start]), size);
    }
    else {
        end -= byteEnd * 4;
                                                                                
        memcpy(sequence, &(bases[*bytes][start]), 4 - start);
        bytes++;
        sequence += (4 - start);
        for (i = byteStart+1; i < byteEnd-1; i++, bytes++, sequence += 4)
            memcpy(sequence, bases[*bytes], 4);
        memcpy(sequence, bases[*bytes], end + 4);
        bytes++;
        bytes -= byteSize;
    }
    return 0;
}

static void
applyNs(char sequence[], Py_ssize_t start, Py_ssize_t end, Py_buffer *nBlocks)
{
    const Py_ssize_t nBlockCount = nBlocks->shape[0];
    const uint32_t* const nBlockPositions = nBlocks->buf;

    Py_ssize_t i;
    for (i = 0; i < nBlockCount; i++) {
        Py_ssize_t nBlockStart = nBlockPositions[2*i];
        Py_ssize_t nBlockEnd = nBlockPositions[2*i+1];
        if (nBlockEnd < start) continue;
        if (end < nBlockStart) break;
        if (nBlockStart < start) nBlockStart = start;
        if (end < nBlockEnd) nBlockEnd = end;
        memset(sequence + nBlockStart - start, 'N', nBlockEnd - nBlockStart);
    }
}

static void
applyMask(char sequence[], Py_ssize_t start, Py_ssize_t end,
          Py_buffer* maskBlocks)
{
    const Py_ssize_t maskBlockCount = maskBlocks->shape[0];
    const uint32_t* const maskBlockPositions = maskBlocks->buf;
    const char diff = 'a' - 'A';

    Py_ssize_t i;
    for (i = 0; i < maskBlockCount; i++) {
        Py_ssize_t j;
        Py_ssize_t maskBlockStart = maskBlockPositions[2*i];
        Py_ssize_t maskBlockEnd = maskBlockPositions[2*i+1];
        if (maskBlockEnd < start) continue;
        if (end < maskBlockStart) break;
        if (maskBlockStart < start) maskBlockStart = start;
        if (end < maskBlockEnd) maskBlockEnd = end;
        for (j = maskBlockStart - start; j < maskBlockEnd - start; j++)
            sequence[j] += diff;
    }
}

static int
blocks_converter(PyObject* object, void* pointer)
{
    const int flag = PyBUF_ND | PyBUF_FORMAT;
    Py_buffer *view = pointer;

    if (object == NULL) goto exit;

    if (PyObject_GetBuffer(object, view, flag) == -1) {
        PyErr_SetString(PyExc_RuntimeError, "blocks have unexpected format.");
        return 0;
    }

    if (view->itemsize != sizeof(uint32_t)
     || (strcmp(view->format, "I") != 0 && strcmp(view->format, "L") != 0 )) {
        PyErr_Format(PyExc_RuntimeError,
                     "blocks have incorrect data type (itemsize %zd, format %s)",
                     view->itemsize, view->format);
        goto exit;
    }
    if (view->ndim != 2) {
        PyErr_Format(PyExc_RuntimeError,
                     "blocks have incorrect rank %d (expected 2)", view->ndim);
        goto exit;
    }
    if (view->shape[1] != 2) {
        PyErr_Format(PyExc_RuntimeError,
                     "blocks should have two columns (found %zd)",
                     view->shape[1]);
        goto exit;
    }
    return Py_CLEANUP_SUPPORTED;

exit:
    PyBuffer_Release(view);
    return 0;
}

static char TwoBit_convert__doc__[] = "convert twoBit data to the DNA sequence, apply blocks of N's (representing unknown sequences) and masked (lower case) blocks, and return the sequence as a bytes object";

static PyObject*
TwoBit_convert(PyObject* self, PyObject* args, PyObject* keywords)
{
    const unsigned char *data;
    Py_ssize_t start;
    Py_ssize_t end;
    Py_ssize_t step;
    Py_ssize_t size;
    Py_ssize_t length;
    Py_buffer nBlocks;
    Py_buffer maskBlocks;
    PyObject *object;
    char *sequence;

    static char* kwlist[] = {"data", "start", "end", "step",
                             "nBlocks", "maskBlocks", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywords, "y#nnnO&O&", kwlist,
                                     &data, &length, &start, &end, &step,
                                     &blocks_converter, &nBlocks,
                                     &blocks_converter, &maskBlocks))
        return NULL;

    size = (end - start) / step;
    object = PyBytes_FromStringAndSize(NULL, size);
    if (!object) goto exit;

    sequence = PyBytes_AS_STRING(object);

    if (step == 1) {
        if (extract(data, length, start, end, sequence) < 0) {
            Py_DECREF(object);
            object = NULL;
            goto exit;
        }
        applyNs(sequence, start, end, &nBlocks);
        applyMask(sequence, start, end, &maskBlocks);
    }
    else {
        Py_ssize_t current, i;
        Py_ssize_t full_start, full_end;
        char* full_sequence;
        if (start <= end) {
            full_start = start;
            full_end = end;
            current = 0;                                 
        }
        else {
            full_start = end + 1;
            full_end = start + 1;
            current = start - end - 1;                                
        }
        full_sequence = PyMem_Malloc((full_end-full_start+1)*sizeof(char));
        full_sequence[full_end-full_start] = '\0';
        if (!full_sequence) {
            Py_DECREF(object);
            object = NULL;
            goto exit;
        }
        if (extract(data, length, full_start, full_end, full_sequence) < 0) {
            PyMem_Free(full_sequence);
            Py_DECREF(object);
            object = NULL;
            goto exit;
        }
        applyNs(full_sequence, full_start, full_end, &nBlocks);
        applyMask(full_sequence, full_start, full_end, &maskBlocks);
        for (i = 0; i < size; current += step, i++)
            sequence[i] = full_sequence[current];
        PyMem_Free(full_sequence);
    }

exit:
    blocks_converter(NULL, &nBlocks);
    blocks_converter(NULL, &maskBlocks);
    return object;
}

static struct PyMethodDef _twoBitIO_methods[] = {
    {"convert",
     (PyCFunction)TwoBit_convert,
     METH_VARARGS | METH_KEYWORDS,
     TwoBit_convert__doc__
    },
    {NULL, NULL, 0, NULL}               
};


static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_twoBitIO",
    "Parser for DNA sequence data in 2bit format",
    -1,
    _twoBitIO_methods,
    NULL,
    NULL,
    NULL,
    NULL
};

PyObject *
PyInit__twoBitIO(void)
{
    return PyModule_Create(&moduledef);
}
