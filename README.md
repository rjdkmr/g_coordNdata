##g_coordNdata
***

###About

This program can be used to extract coordinates from the GROMACS MD 
Trajectory. The output format of coordinate data file is as follows:
<pre><code>--------------------------------------------------------------------
[number of frames]       [number of atoms]                          
X11 Y11 Z11 X12 Y12 Z12 X13 Y13 Z13 ... ... ... ... ... X1n Y1n Z1n 
X21 Y21 Z21 X22 Y22 Z22 X23 Y23 Z23 ... ... ... ... ... X2n Y2n Z2n 
X31 Y31 Z31 X32 Y32 Z32 X33 Y33 Z33 ... ... ... ... ... X3n Y3n Z3n 
... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... ... 
Xm1 Ym1 Zm1 Xm2 Ym2 Zm2 Xm3 Ym3 Zm3 ... ... .. ... ... Xmn Ymn Zmn  
--------------------------------------------------------------------
</pre></code>
where "m" is number of frames and "n" is number of atoms. This file could be
used in other program or scripts for MATLAB, R or Python for further analysis
and coordinates transformation.                 

The program could read coordinates data from an external file with the above
described file format to create a GROMACS format MD trajectory. This is
useful when after analysis and coordinates transformations visualization of
trajectory is required.

***

###Requirements
To compile and install, GROMACS libraries <code> libgmx, libmd, libgmxana </code> are required.
***

###Download
<pre><code>git clone https://github.com/rjdkmr/g_coordNdata
</code></pre>
***

###Installation
<pre><code>cd g_coordNdata
mkdir build
cd build
cmake ..  -DGMX_PATH=/opt/gromacs
make
make install
</code></pre>

Directory <code>/opt/gromacs</code> should contains <code>include</code> and <code> lib </code> directories. If these directories are in seprate locations, use followings:
<pre><code>cmake ..  -DGMX_LIB=/path/to/lib -GMX_INCLUDE=/path/to/include
</code></pre>

If fftw library <code> libfftw3f.so or libfftw3f.a </code> are not present in standard locations:
<pre><code>-DFFTW_LIB=/path/to/fftw3/lib</code></pre>
***

###Usage
<pre><code>g_coordNdata -h
</code></pre>
***
