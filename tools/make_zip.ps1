# Build the zip that can be sent to someone with a private repo.
# Run from the repository root:  pwsh tools/make_zip.ps1

$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $PSScriptRoot
$stage = Join-Path $env:TEMP "chemdraft-stage"
$zip = Join-Path ([Environment]::GetFolderPath("Desktop")) "chemdraft.zip"

Remove-Item $stage -Recurse -Force -ErrorAction SilentlyContinue
New-Item -ItemType Directory -Force $stage | Out-Null

$include = @(
    "chemdraft", "skills", "examples", "tools", "tests",
    "pyproject.toml", "README.md", "VALIDATION.md", "SETUP.md",
    "install.cmd", "update.cmd", ".gitignore", ".gitattributes"
)

foreach ($item in $include) {
    $source = Join-Path $root $item
    if (Test-Path $source) { Copy-Item -Recurse -Force $source $stage }
}

Get-ChildItem $stage -Recurse -Directory |
    Where-Object { $_.Name -in @("__pycache__", ".pytest_cache") } |
    Remove-Item -Recurse -Force -ErrorAction SilentlyContinue

if (Test-Path $zip) { Remove-Item $zip -Force }
Compress-Archive -Path (Join-Path $stage "*") -DestinationPath $zip -CompressionLevel Optimal

$size = [math]::Round((Get-Item $zip).Length / 1MB, 2)
$files = (Get-ChildItem $stage -Recurse -File).Count
Write-Host "built $zip  ($size MB, $files files)"
