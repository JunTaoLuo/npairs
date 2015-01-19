package timing;

import java.lang.management.*;


public class Time {
	
/** get CPU time in nanoseconds. */
public long getCpuTime() {
	ThreadMXBean bean = ManagementFactory.getThreadMXBean();
	return bean.isCurrentThreadCpuTimeSupported() ?
			bean.getCurrentThreadCpuTime() : 0L;
	
}

/** get CPU time in seconds. */
public double getCpuTimeSec() {
	return getCpuTime() / Math.pow(10, 9);		
}

/** Get user time in nanoseconds. */
public long getUserTime( ) {
    ThreadMXBean bean = ManagementFactory.getThreadMXBean( );
    return bean.isCurrentThreadCpuTimeSupported( ) ?
        bean.getCurrentThreadUserTime( ) : 0L;
}

/** Get user time in seconds. */
public double getUserTimeSec() {
	return getUserTime() / Math.pow(10, 9);
}

/** Get system time in nanoseconds. */
public long getSystemTime( ) {
    ThreadMXBean bean = ManagementFactory.getThreadMXBean( );
    return bean.isCurrentThreadCpuTimeSupported( ) ?
        (bean.getCurrentThreadCpuTime( ) - bean.getCurrentThreadUserTime( )) : 0L;
}

/** Get system time in seconds. */
public double getSystemTimeSec() {
	return getSystemTime() / Math.pow(10, 9);
}

/** Get JVM CPU time in milliseconds (using undocumented internal com.sun.OperatingSystemMXBean) */
public long getJVMCpuTime( ) {
    OperatingSystemMXBean bean =
        ManagementFactory.getOperatingSystemMXBean( );
    if ( ! (bean instanceof
        com.sun.management.OperatingSystemMXBean) )
        return 0L;
    return ((com.sun.management.OperatingSystemMXBean)bean)
        .getProcessCpuTime( );
}

}
