$root = (Resolve-Path 'results/sub6_planA/stage1').Path
$watch = [System.IO.FileSystemWatcher]::new($root, 'ANCHOR_*.json')
$watch.IncludeSubdirectories = $true
$watch.NotifyFilter = [System.IO.NotifyFilters]::LastWrite
$watch.EnableRaisingEvents = $true
$subscription = Register-ObjectEvent -InputObject $watch -EventName Changed -SourceIdentifier Stage1ResultChanged
try {
    while ($true) {
        $rows = @(Get-ChildItem -Path $root -Directory -Filter 'ANCHOR_*' | ForEach-Object {
            $p = Join-Path $_.FullName ($_.Name + '.json')
            if (Test-Path -LiteralPath $p) {
                try { Get-Content -LiteralPath $p -Raw | ConvertFrom-Json } catch { $null }
            }
        } | Where-Object { $_ })
        $failed = @($rows | Where-Object { $_.status -eq 'FAILED' })
        $done = @($rows | Where-Object { $_.status -eq 'SOLVED_EXPORTED' })
        if ($failed.Count -gt 0) {
            Write-Output "FAILED: $($failed.run_id -join ', ')"
            break
        }
        if ($done.Count -eq 10) {
            Write-Output 'ALL 10 SOLVED_EXPORTED'
            break
        }
        $event = Wait-Event -SourceIdentifier Stage1ResultChanged -Timeout 1800
        if ($null -ne $event) { Remove-Event -EventIdentifier $event.EventIdentifier }
    }
} finally {
    Unregister-Event -SourceIdentifier Stage1ResultChanged
    $watch.Dispose()
}
