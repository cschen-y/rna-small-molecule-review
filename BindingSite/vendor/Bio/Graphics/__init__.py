







"""Bio.Graphics offers several graphical outputs, all using ReportLab."""


try:
    import reportlab as r

    del r
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Please install ReportLab if you want "
        "to use Bio.Graphics. You can find ReportLab at "
        "http://www.reportlab.com/software/opensource/"
    ) from None







def _write(drawing, output_file, format, dpi=72):
    """Standardize output to files (PRIVATE).

    Writes the provided drawing out to a file in a prescribed format.

      - drawing - suitable ReportLab drawing object.
      - output_file - a handle to write to, or a filename to write to.
      - format - String indicating output format, one of PS, PDF, SVG,
        or provided the ReportLab renderPM module is installed,
        one of the bitmap formats JPG, BMP, GIF, PNG, TIFF or TIFF.
        The format can be given in any case.
      - dpi - Resolution (dots per inch) for bitmap formats.

    No return value.
    """
    from reportlab.graphics import renderPDF
    from reportlab.graphics import renderPS
    from reportlab.graphics import renderSVG

    try:
        from reportlab.graphics import renderPM
    except ImportError:
        
        
        
        renderPM = None

    formatdict = {
        "PS": renderPS,
        "EPS": renderPS,
        
        
        "PDF": renderPDF,
        "SVG": renderSVG,
        "JPG": renderPM,
        "BMP": renderPM,
        "GIF": renderPM,
        "PNG": renderPM,
        "TIFF": renderPM,
        "TIF": renderPM,
    }
    try:
        
        
        drawmethod = formatdict[format.upper()]  
    except (KeyError, AttributeError):
        raise ValueError(
            f"Output format should be one of {', '.join(formatdict)}"
        ) from None

    if drawmethod is None:
        
        
        from Bio import MissingPythonDependencyError

        raise MissingPythonDependencyError("Please install ReportLab's renderPM module")

    if drawmethod == renderPM:
        
        return drawmethod.drawToFile(drawing, output_file, format, dpi=dpi)
    else:
        return drawmethod.drawToFile(drawing, output_file)
