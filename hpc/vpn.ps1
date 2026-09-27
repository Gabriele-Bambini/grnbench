# Bocconi VPN (Palo Alto GlobalProtect) on Windows, for people and agents.
#
# The HPC login node resolves only inside the VPN. GlobalProtect has no Windows
# CLI; the scriptable action is restarting its PanGPS service, which rebuilds the
# tunnel while the portal still holds a cached sign-in cookie. When the cookie has
# expired, GlobalProtect opens its SSO window and the PERSON signs in: no script
# does that part.
#
# PanGPS runs as SYSTEM. A one-time grant (admin shell) lets your own account
# query, start and stop that one service, so reconnect needs no elevation after.
#
#   powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 status      # VPN up / VPN down
#   powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 reconnect   # no-op if up; else restart PanGPS, wait up to 90 s
#   powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 setup       # print the one-time admin command
#   (admin) ... -File hpc/vpn.ps1 grant <SID>                                   # what `setup` prints; idempotent
#
# Exit codes (reconnect): 0 up and the login node resolves, 1 not so after 90 s (sign in once in the
# GlobalProtect window), 2 the one-time rights grant is missing.

param(
    [Parameter(Position = 0)][ValidateSet("status", "reconnect", "setup", "grant")][string]$cmd = "status",
    [Parameter(Position = 1)][string]$sid = "",
    [string]$LoginHost = "slogin.hpc.unibocconi.it"
)

function Get-VpnAdapter {
    Get-NetAdapter -ErrorAction SilentlyContinue |
        Where-Object { $_.InterfaceDescription -match "PANGP" }
}

function Test-LoginHost {
    # The adapter can come up seconds before the login node resolves (measured 2026-09-27):
    # "up" only counts once the name resolves.
    try { [void][System.Net.Dns]::GetHostAddresses($LoginHost); return $true } catch { return $false }
}

function Show-Status {
    $ad = Get-VpnAdapter
    if ($ad -and $ad.Status -eq "Up") {
        $res = if (Test-LoginHost) { "resolves" } else { "does not resolve yet" }
        Write-Output "VPN up ($($ad.Name)); $LoginHost $res"
        return
    }
    Write-Output "VPN down"
}

if ($cmd -eq "status") { Show-Status; exit 0 }

if ($cmd -eq "setup") {
    # Run NON-elevated: the SID must be the account that will run reconnect, which
    # can differ from the admin account that applies the grant.
    $me = ([System.Security.Principal.WindowsIdentity]::GetCurrent()).User.Value
    $script = (Resolve-Path $PSCommandPath).Path
    Write-Output "One-time setup. Open PowerShell as administrator and run:"
    Write-Output ""
    Write-Output "  powershell -NoProfile -ExecutionPolicy Bypass -File `"$script`" grant $me"
    Write-Output ""
    Write-Output "It adds one entry to the PanGPS service's permissions: query, list dependents, start and stop, for"
    Write-Output "this account ($me) only. Then run: powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 reconnect"
    exit 0
}

if ($cmd -eq "grant") {
    # Local/AD accounts are S-1-5-21-..., Microsoft Entra ID accounts S-1-12-1-...
    if ($sid -notmatch '^S-1-(5-21|12-1)(-\d+){4}$') {
        Write-Output "grant needs your account SID (S-1-5-21-... or S-1-12-1-...). Run 'hpc/vpn.ps1 setup' non-elevated: it prints the exact command."
        exit 2
    }
    $admin = ([Security.Principal.WindowsPrincipal][Security.Principal.WindowsIdentity]::GetCurrent()).IsInRole(
        [Security.Principal.WindowsBuiltInRole]::Administrator)
    if (-not $admin) { Write-Output "grant must run in an administrator PowerShell."; exit 2 }
    # LC query status, SW enumerate dependents (Restart-Service reads them), RP start, WP stop.
    $ace = "(A;;LCSWRPWP;;;$sid)"
    $raw = ((sc.exe sdshow PanGPS) -join '').Trim()
    if (-not $raw.StartsWith("D:")) { Write-Output "cannot read the PanGPS security descriptor: $raw"; exit 2 }
    if ($raw.Contains($ace)) { Write-Output "already granted: $ace"; exit 0 }
    # The new ACE goes at the end of the DACL, before the SACL ("S:") if there is one.
    $i = $raw.IndexOf("S:")
    if ($i -lt 0) { $new = $raw + $ace } else { $new = $raw.Substring(0, $i) + $ace + $raw.Substring($i) }
    sc.exe sdset PanGPS $new
    if ($LASTEXITCODE -ne 0) { exit 2 }
    Write-Output "granted $ace on PanGPS. A GlobalProtect update may reset it: run setup again if reconnect says access denied."
    exit 0
}

# reconnect. A working tunnel is never torn down: restarting PanGPS cuts every
# session using it (other agents' ssh, a running rsync). Restart only when the
# adapter is down or the login node does not resolve.
$ad = Get-VpnAdapter
if ($ad -and $ad.Status -eq "Up" -and (Test-LoginHost)) {
    Write-Output "VPN already up ($($ad.Name)); $LoginHost resolves. Nothing restarted."
    exit 0
}
try { Restart-Service PanGPS -ErrorAction Stop }
catch {
    Write-Output "FAILED: this account cannot restart the GlobalProtect service (PanGPS): $($_.Exception.Message)"
    Write-Output "Most often the one-time grant is missing or was reset by a GlobalProtect update."
    Write-Output "One-time setup: powershell -NoProfile -ExecutionPolicy Bypass -File hpc/vpn.ps1 setup"
    exit 2
}

# Wait up to 90 s. An expired portal cookie shows up here as a timeout while
# GlobalProtect waits for the person in its sign-in window.
$upAt = $null
foreach ($i in 1..45) {
    Start-Sleep -Seconds 2
    $ad = Get-VpnAdapter
    if (-not ($ad -and $ad.Status -eq "Up")) { continue }
    if ($null -eq $upAt) { $upAt = $i * 2 }
    if (Test-LoginHost) { Write-Output "VPN up after $($upAt)s; $LoginHost resolves after $($i * 2)s"; exit 0 }
}
if ($null -ne $upAt) {
    Write-Output "VPN adapter up after $($upAt)s but $LoginHost still does not resolve after 90s."
} else {
    Write-Output "VPN still down after 90s: the portal wants an interactive sign-in. Open GlobalProtect and connect once (portal vpnstudents.unibocconi.it); later reconnects reuse that sign-in."
}
exit 1
