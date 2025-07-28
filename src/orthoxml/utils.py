# utils.py

import bz2
import gzip
import os
from io import BytesIO
from importlib import resources
from lxml import etree
from .logger import logger

# File opening. This is based on the example on SO here:
# http://stackoverflow.com/a/26986344
fmagic = {b'\x1f\x8b\x08': gzip.open,
          b'\x42\x5a\x68': bz2.BZ2File}


def auto_open(fn, *args, **kwargs):
    """function to open regular or compressed files for read / write.

    This function opens files based on their "magic bytes". Supports bz2
    and gzip. If it finds neither of these, presumption is it is a
    standard, uncompressed file.

    Example::

        with auto_open("/path/to/file/maybe/compressed", mode="rb") as fh:
            fh.read()

        with auto_open("/tmp/test.txt.gz", mode="wb") as fh:
            fh.write("my big testfile")

    :param fn: either a string of an existing or new file path, or
        a BytesIO handle
    :param **kwargs: additional arguments that are understood by the
        underlying open handler
    :returns: a file handler
    """
    if isinstance(fn, BytesIO):
        return fn

    if os.path.isfile(fn) and os.stat(fn).st_size > 0:
        with open(fn, 'rb') as fp:
            fs = fp.read(max([len(x) for x in fmagic]))
        for (magic, _open) in fmagic.items():
            if fs.startswith(magic):
                return _open(fn, *args, **kwargs)
    else:
        if fn.endswith('gz'):
            return gzip.open(fn, *args, **kwargs)
        elif fn.endswith('bz2'):
            return bz2.BZ2File(fn, *args, **kwargs)

    return open(fn, *args, **kwargs)

def validate_xml(xml_file_path: str, orthoxml_version: str) -> bool:
    """
    Stream-validate an OrthoXML file against the specified schema version.

    :param xml_file_path: Path to the OrthoXML file.
    :param orthoxml_version: The OrthoXML version as a string.
    :return: True if valid, False otherwise.
    """
    try:
        # Load XSD schema from package resources
        schema_filename = f'orthoxml-{orthoxml_version}.xsd'
        with resources.files('orthoxml.schemas').joinpath(schema_filename).open('rb') as schema_file:
            schema_root = etree.XML(schema_file.read())
            schema = etree.XMLSchema(schema_root)

        # Stream-parse and validate the file
        for _, elem in etree.iterparse(xml_file_path, schema=schema, recover=False):
            elem.clear()

        return True

    except etree.XMLSyntaxError as e:
        logger.error(f"Validation failed: XML file '{xml_file_path}' is not valid for schema {schema_filename}: {e}")
        return False

    except etree.XMLSchemaParseError as e:
        logger.error(f"Invalid XML Schema '{schema_filename}': {e}")
        raise

    except Exception as e:
        logger.exception(f"Unexpected error during validation of '{xml_file_path}': {e}")
        raise