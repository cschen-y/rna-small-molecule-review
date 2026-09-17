$vendor = Join-Path $PSScriptRoot "vendor"
if ($env:PYTHONPATH) {
    $env:PYTHONPATH = "$vendor;$env:PYTHONPATH"
} else {
    $env:PYTHONPATH = $vendor
}
& "D:\code_tool\anaconda\python.exe" @args
exit $LASTEXITCODE
