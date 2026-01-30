"""
Resource detection utilities for parallel processing.

Provides cgroups-aware detection of CPU and memory resources for
use in parallel sample processing.
"""

import os
import logging
from typing import Tuple, Optional

logger = logging.getLogger(__name__)


def detect_available_cpus() -> int:
    """
    Detect available CPUs, respecting cgroups/SLURM allocation.
    
    Returns:
        Number of available CPU cores.
    """
    try:
        # sched_getaffinity respects cgroups/SLURM task allocation
        return len(os.sched_getaffinity(0))
    except AttributeError:
        # macOS doesn't have sched_getaffinity
        pass
    
    # Fallback to os.cpu_count()
    return os.cpu_count() or 4


def _get_cgroup_memory_limit_gb() -> Optional[float]:
    """
    Read memory limit from cgroups (SLURM/Docker).
    
    Returns:
        Memory limit in GB, or None if not in a cgroup.
    """
    # Try cgroups v2 first
    try:
        with open('/sys/fs/cgroup/memory.max', 'r') as f:
            val = f.read().strip()
            if val != 'max':
                return int(val) / (1024 ** 3)
    except (FileNotFoundError, PermissionError, ValueError):
        pass
    
    # Try cgroups v1
    try:
        with open('/sys/fs/cgroup/memory/memory.limit_in_bytes', 'r') as f:
            val = int(f.read().strip())
            # Ignore if set to max (huge value close to max int)
            if val < 10 ** 18:  # Less than ~1 exabyte
                return val / (1024 ** 3)
    except (FileNotFoundError, PermissionError, ValueError):
        pass
    
    return None


def _get_meminfo_available_gb() -> Optional[float]:
    """
    Read available memory from /proc/meminfo (Linux only).
    
    Returns:
        Available memory in GB, or None on non-Linux.
    """
    try:
        with open('/proc/meminfo', 'r') as f:
            for line in f:
                if line.startswith('MemAvailable:'):
                    kb = int(line.split()[1])
                    return kb / (1024 ** 2)  # KB to GB
    except (FileNotFoundError, PermissionError, ValueError):
        pass
    
    return None


def detect_available_memory_gb() -> float:
    """
    Detect available memory, respecting cgroups limits.
    
    Checks in order:
    1. cgroups memory limit (SLURM/Docker)
    2. /proc/meminfo MemAvailable (Linux)
    3. psutil (cross-platform, if installed)
    4. Fallback to 32 GB
    
    Returns:
        Available memory in GB.
    """
    # Check cgroups first (SLURM/Docker environments)
    cgroup_limit = _get_cgroup_memory_limit_gb()
    if cgroup_limit is not None:
        logger.debug(f"Detected cgroup memory limit: {cgroup_limit:.1f} GB")
        return cgroup_limit
    
    # Check /proc/meminfo (Linux)
    meminfo = _get_meminfo_available_gb()
    if meminfo is not None:
        logger.debug(f"Detected available memory from meminfo: {meminfo:.1f} GB")
        return meminfo
    
    # Try psutil as last resort
    try:
        import psutil
        available = psutil.virtual_memory().available / (1024 ** 3)
        logger.debug(f"Detected available memory from psutil: {available:.1f} GB")
        return available
    except ImportError:
        pass
    
    # Conservative fallback
    logger.debug("Could not detect memory, using fallback: 32 GB")
    return 32.0


def detect_resources() -> Tuple[int, float]:
    """
    Detect available CPU and memory resources.
    
    Returns:
        Tuple of (cpu_count, memory_gb).
    """
    cpus = detect_available_cpus()
    memory_gb = detect_available_memory_gb()
    return cpus, memory_gb


def calculate_auto_parallel_samples(
    threads: int,
    memory_gb: float,
    memory_per_sample: float = 2.0,
    min_threads_per_sample: int = 4,
) -> int:
    """
    Calculate optimal number of parallel samples based on resources.
    
    Args:
        threads: Total available threads.
        memory_gb: Available memory in GB.
        memory_per_sample: Expected peak memory per sample in GB.
        min_threads_per_sample: Minimum threads per sample for efficiency.
    
    Returns:
        Recommended number of parallel samples.
    """
    # Memory constraint: don't exceed available RAM
    mem_limit = int(memory_gb / memory_per_sample)
    
    # CPU constraint: ensure each sample has enough threads for rayon
    cpu_limit = threads // max(1, min_threads_per_sample)
    
    # Take the more restrictive limit
    parallel = max(1, min(mem_limit, cpu_limit))
    
    logger.debug(f"Auto parallel samples: mem_limit={mem_limit}, cpu_limit={cpu_limit}, result={parallel}")
    return parallel
