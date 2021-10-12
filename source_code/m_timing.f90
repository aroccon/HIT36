module m_timing
integer*8 :: cpu0, cpu1, cpu2, dcpu
integer*4 :: cpu_sec, cpu_min, cpu_hrs, cpu_min_total
integer*4 :: job_runlimit=9999


contains



subroutine m_timing_init

call system_clock(cpu0,dcpu)

return
end subroutine m_timing_init







subroutine m_timing_check
use m_openmpi
implicit none

call system_clock(cpu1,dcpu)
cpu_sec = (cpu1-cpu0)/dcpu
cpu_min = cpu_sec/60;   cpu_min_total = cpu_min
cpu_hrs = cpu_min/60
cpu_min = mod(cpu_min,60)
cpu_sec = mod(cpu_sec,60)

return
end subroutine m_timing_check







end module m_timing
