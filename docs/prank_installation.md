Installation of PRANK
=====================

*   [Pre-compiled executables for Windows](#pre-compiled-executables-for-windows)
    *   [Hints for using PRANK on Windows](#hints-for-using-prank-on-windows)
*   [Pre-compiled executables for Mac OSX](#pre-compiled-executables-for-mac-osx)
*   [Pre-compiled executables for Linux](#pre-compiled-executables-for-linux)
*   [Using PRANK with Docker](#using-prank-with-docker)
*   [Installation of PRANK from source code](#installation-of-prank-from-source-code)
    *   [Installation of helper tools](#installation-of-helper-tools)

   
[Back to PRANK home.](../README.md)  
 

### Pre-compiled executables for Windows

A pre-compiled PRANK executable for Windows is available in the [binaries folder](../binaries/). The executable runs within the [Cygwin environment](https://www.cygwin.com/). Please follow [the instructions](https://cygwin.com/install.html) and install the default (or a more complete) setup of Cygwin.

Assuming that Cygwin is installed at “C:\\cygwin64”, it is recommended to drop the file prank.cygwin.\*.zip file to the directory “C:\\cygwin64\\home\\username”. You can either unpack the zip file using a graphicla program, or (if you have installed Cygwin’s program “unzip”) open a Cygwin64 Terminal and type the commands:  
 
```
cd
unzip prank.cygwin.*.zip
```
   
You can then start PRANK within the Cygwin terminal using the command:

```
./prank.cygwin/prank/bin/prank
```

The executable has been compiled and tested on Windows 10 (64bit) system.

#### Hints for using PRANK on Windows

PRANK is used through commands given in a Cygwin terminal. For those not familiar with moving around files, the following workflow may be useful.

*   open “Cygwin64 Terminal”
*   in the terminal, type commands

```
mkdir analysis
cd analysis
ln -s ~/prank.cygwin/prank/bin/prank.exe .
explorer .
```

*   drag and drop your alignment input file(s) to the Explorer window that was opened
*   in a Cygwin terminal, run PRANK as ./prank -d=myfile \[other options\]
*   once finished, the result files will appear in the Explorer window

### Pre-compiled executables for Mac OSX

The easiest way to install PRANK (together with its helper programs MAFFT and Exonerate) is to use the package manager [Homebrew](https://brew.sh): 

```
brew install brewsci/bio/prank
```

Alternatively, you can use a pre-compiled package that is available in the [binaries folder](../binaries/). The zip archive contains standalone executables for PRANK and the helper tools BppAncestor, Exonerate and MAFFT. It has been tested with MacOS version 12.1.

The package file can be downloaded and unpacked using graphical software, or on the command line:

```
curl -O http://wasabiapp.org/download/prank/prank.osx64.170427.zip
unzip prank.osx64.170427.zip
```

Then, remove the MacOS quarantine flag from the downloaded files (prevents “prank cannot be opened” error message) and test the executable:

```
sudo xattr -r -d com.apple.quarantine prank-msa
cd prank-msa
./prank -d=../my_test_file.fa
```

For the ease of use, it is recommended to add the prank-msa directory to the system path or copy its content to a path directory (such as ~/bin or /usr/local/bin). This can be done with commands similar to these:

```
echo "export PATH=$PATH:/Users/$USER/Downloads/prank-msa" >> ~/.bash_profile
```

or 

```
cp -R prank-msa/* /usr/local/bin
```

Note that if the executables are moved or copied, the directory structure of the helper tools (e.g. mafftdir/bin) has to be retained.

If you want to use a separate, system wide installation for any of the PRANK helper tools (instead of the version bundled with the package), you can just remove it from the prank directory (e.g. rm -rf mafft\*  or rm exonerate). Installing the helper tools globally is described in a dedicated section below.

### Pre-compiled executables for Linux

The pre-compiled executables are available in the [binaries folder](../binaries/).

The PRANK executable for Linux has been compiled on Red Hat 5.10 (64bit) and confirmed to work also on Ubuntu 13.10 and CentOS 6.5. The package contains the PRANK executable and its manual page (read it with command man ./prank/prank.1) as well as BppAncestor, Exonerate and MAFFT executables and the necessary files for those to work correctly.

Using the command line, the PRANK executable for Linux can be downloaded and unpacked with the commands:

```
mkdir ~/programs
cd ~/programs
wget http://wasabiapp.org/download/prank/prank.linux64.140603.tgz
tar prank.linux64.140603.tgz
./prank/bin/prank
```

For the ease of use, it is recommended to add the directory prank/bin to the system path or copy its content to a such directory (such as ~/bin or /usr/local/bin). This can be done with commands similar to these:

```
echo "export PATH=$PATH:/home/$USER/programs/prank/bin" >> ~/.bash_profile
```

or

```
cp -R /home/$USER/programs/prank/bin/* ~/bin/
```

Note that if the executables of the bundled version are moved or copied, the directory structure (bin/lib) and the location of library dependencies have to be retained.

### Using PRANK with Docker

PRANK for Linux is provided also as Docker image that can be used on all computer platforms supporting [Docker](https://www.docker.com/).  
 

**Download PRANK image:**

```
docker pull ariloytynoja/prank
docker tag ariloytynoja/prank prank
```

**Run the image:**

```
docker run --rm -v $path\_to\_current\_directory:/data prank -d=infile.fas +F
```

**Alternative: Create a helper script (for Linux or OSX) and run that:**

```
cat > prank.sh << EOF 
#!/bin/bash
docker run --rm -v \`pwd\`:/data prank "\\$@"
EOF

chmod +x prank.sh 

./prank.sh -d=infile.fas +F
```

### Installation of PRANK from source code

You can download a snapshot of the PRANK code from the [binaries folder](../binaries/) and get the latest version using git as explained on the [source code page](https://github.com/ariloytynoja/prank-msa). On command line, the process looks like this (tested on Linux and MacOSX):

```
wget http://wasabiapp.org/download/prank/prank.source.140603.tgz
```
or
```
curl -O http://wasabiapp.org/download/prank/prank.source.140603.tgz
```

```
tar xvzf prank.source.140603.tgz
cd prank-msa/src
make
./prank

[ cp prank ~/bin ]
[ sudo cp prank /usr/bin ]
```

or

```
git clone https://github.com/ariloytynoja/prank-msa.git 
cd prank-msa/src
make
./prank

[ cp prank ~/bin ]
[ sudo cp prank /usr/bin ]
```

On Windows, the PRANK source code has been verified to compile on the [Cygwin environment](http://cygwin.com/).

#### Installation of helper tools

For the fast guide tree estimation, alignment anchoring and ancestral reconstruction to work, you need to have MAFFT, Exonerate and BppAncestor installed and on the execution path. On many popular Linux distributions MAFFT and Exonerate are available in the software repository and can be installed pre-compiled. On Ubuntu, that is done with the following commands:

```
sudo apt-get install mafft
sudo apt-get install exonerate
```

On MacOS, MAFFT and Exonerate can be installed with a package management system [Homebrew](https://brew.sh):

```
brew install mafft
brew install brewsci/bio/exonerate
```

Alternatively, you can download the MAFFT source code from [http://mafft.cbrc.jp/alignment/software](http://mafft.cbrc.jp/alignment/software), the Exonerate source code from [https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) and follow the instructions for their installation.

BppSuite is provided as pre-compiled binaries for several platforms at [https://github.com/BioPP/bppsuite](https://github.com/BioPP/bppsuite) and is also available in repositories of some popular Linux distros. Unfortunately, the version of BppAncestor provided there probably contains a critical bug and cannot be used. The commands below compile the program suite against the latest version of the Bpp libraries.

```
mkdir bppsuite
cd bppsuite/
bpp_dir=`pwd`

mkdir sources 
cd sources/

git clone http://biopp.univ-montp2.fr/git/bpp-core
git clone http://biopp.univ-montp2.fr/git/bpp-seq
git clone http://biopp.univ-montp2.fr/git/bpp-phyl
git clone http://biopp.univ-montp2.fr/git/bppsuite


cd bpp-core/
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_LIBRARY_PATH=$bpp_dir/lib -DCMAKE_INCLUDE_PATH=$bpp_dir/include -D BUILD_TESTING=FALSE -DCMAKE_EXE_LINKER_FLAGS=-static ./
make
make install

cd ../bpp-seq/
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_LIBRARY_PATH=$bpp_dir/lib -DCMAKE_INCLUDE_PATH=$bpp_dir/include -D BUILD_TESTING=FALSE -DCMAKE_EXE_LINKER_FLAGS=-static ./
make
make install

cd ../bpp-phyl/
cmake -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_LIBRARY_PATH=$bpp_dir/lib -DCMAKE_INCLUDE_PATH=$bpp_dir/include -D BUILD_TESTING=FALSE -DCMAKE_EXE_LINKER_FLAGS=-static ./
make
make install

cd ../bppsuite/
cmake . -DCMAKE_PREFIX_PATH=$bpp_dir -DCMAKE_INSTALL_PREFIX=$bpp_dir -DCMAKE_EXE_LINKER_FLAGS=-static -DBUILD_STATIC=true
make
make install

cd ../../bin/
ls | xargs strip
```
