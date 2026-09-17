import sys
from distutils.core import Extension, setup


extra_compile_args = [] if sys.platform == "win32" else ["-D_GLIBCXX_USE_CXX11_ABI=0"]
module = Extension(
    "alignment_C",
    sources=["alignment.cpp"],
    extra_compile_args=extra_compile_args,
)

setup(name="alignment_C", version="1.0.0", ext_modules=[module])
