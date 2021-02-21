FFTW
^^^^

(information taken from `<http://micro.stanford.edu/wiki/Install_FFTW3>`_)

If you do not have sudo rights you want to install the libraries in the user space. First create a couple of directories for that purpose::
  
   > mkdir $HOME/usr
   > mkdir $HOME/soft
 
To install FFTW3, download the package from the `FFTW3 download <http://www.fftw.org/download.html>`_ page and decompress it::
  
   > cd ~/soft
   > wget http://www.fftw.org/fftw-3.3.8.tar.gz 
   > tar -zxvf fftw-3.3.8.tar.gz 
   > cd fftw-3.3.8 


**Ubuntu only:** If you want to install FFTW3 (serial version) in your local Ubuntu just run::
  
   > sudo apt-get install libfftw3-dev libfftw3-doc
 
 
However the MPI version will not be available. If you want to have the MPI version follow the instructions in the other sections.

FFTW serial version
^^^^^^^^^^^^^^^^^^^

Configure, make and install::
  
   > ./configure --prefix=$HOME/usr --enable-shared=yes
   > make
   > make install
   


The following files will be installed in::
  
   > ~/usr/include/fftw3.h
   > ~/usr/include/fftw3.f
   > ~/usr/lib/libfftw3.a
   > ~/usr/lib/libfftw3.la
   > ~/usr/lib/libfftw3.so


