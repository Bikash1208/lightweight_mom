"""Time one library task (4 states) at F = 2n via the build script's own solve_point. usage: time_f2n.py MESH THETA PHI"""
import ctypes
import ctypes.wintypes as wt
import sys
import time

sys.path.insert(0, r"C:\lightweight_mom\scripts\research")
sys.path.insert(0, r"C:\lightweight_mom")
import phase20_r4_build_library as b  # noqa: E402


class PMC(ctypes.Structure):
    _fields_ = [("cb", wt.DWORD), ("PageFaultCount", wt.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t),
                ("WorkingSetSize", ctypes.c_size_t), ("a", ctypes.c_size_t), ("b", ctypes.c_size_t),
                ("c", ctypes.c_size_t), ("d", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t),
                ("PeakPagefileUsage", ctypes.c_size_t)]


n, th, ph = int(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
t = time.time()
b._solve_point((10.5, th, ph, n))  # bypasses the mesh-120 point cache on purpose
dt = time.time() - t
p = PMC(); p.cb = ctypes.sizeof(p)
k = ctypes.WinDLL("kernel32"); k.GetCurrentProcess.restype = wt.HANDLE
k.K32GetProcessMemoryInfo.argtypes = [wt.HANDLE, ctypes.POINTER(PMC), wt.DWORD]
k.K32GetProcessMemoryInfo(k.GetCurrentProcess(), ctypes.byref(p), p.cb)
print(f"mesh {n} F={b.floquet_order(n)} theta {th} phi {ph}: {dt:.1f} s, peak {p.PeakWorkingSetSize/1e9:.2f} GB")
