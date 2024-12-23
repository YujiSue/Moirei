# Moirei



# How to install

## Dependency
This applictioion requires the following programs to build and inastall.  

* Python
* CMake
* Build tools
  - For Unix, Linux, Mac  
    GCC
    Clang 
  - For Windows  
    cl.exe (MSVC)

## Run installation command

```sh
git clone https://github.com/YujiSue/Moirei.git
cd Moirei
cmake -DINSTALL_SLIB=ON -S . -B build
cmake --build build
cmake --install build --prefix path-to-app
```

Please replace "path-to-app" to your prefered location.
This application requires some libraries of [slib](https://github.com/YujiSue/slib), and the above command will install the **slib** automatically if need. When you have installed the **slib** already, the option "-DINSTALL_SLIB=ON" will be ignored.

\* If you are restricted to access the system directories (i.e. when you are using a rental server), set "-DBUILD_SLIB=ON" option instead of "-DINSTALL_SLIB=ON". The option allows you to build required libraries in the Moirei's build directory. In that case, do NOT forget to set the environmental variable to indicate the location of **slib** libraries. 


## Another way to install

You can also install this app via my another python library [ysngs](https://github.com/YujiSue/ysngs).

Install the library and call the installation function.

```sh
pip install git+https://github.com/YujiSue/ysngs.git
```
```py
from ysngs import AppManager
ysapps = AppManager()
ysapps.install('moirei')
```

# How to use

You can check all the commands using the help option.
```sh
moirei --help
```

To view the detailed options for each command, set command and run help. 
```sh
moirei command --help
```

## Command list


## Trial 
Trial usage of each function is available in the [Google colab notebook]().

# Preparation of dataset
## Prebuild dataset
The prebuild dataset for human and worm (C. elegans) are available from the following link.

> Human dataset (GRCh38.p14 derived from [RefSeq](https://ftp.ncbi.nlm.nih.gov/genomes/refseq/))  
> [Genomic binary](https://firebasestorage.googleapis.com/v0/b/publicstorage-3ef6a.appspot.com/o/human.bin?alt=media&token=66d4d1f7-57c0-47e9-9b3e-2778a58ef9a1) (Only chr 1-22,X,Y+M; The other short linkage groups and contigs are omitted.)  
> [Annotation database]()  
  

> Nematode worm dataset (C. elegans; WS291 derived from [WormBase](https://downloads.wormbase.org/species/c_elegans/))  
> [Genomic binary](https://firebasestorage.googleapis.com/v0/b/publicstorage-3ef6a.appspot.com/o/worm.bin?alt=media&token=5b3867c1-9a35-46c3-aefd-4a792b697188)  
> [Annotation database](https://firebasestorage.googleapis.com/v0/b/publicstorage-3ef6a.appspot.com/o/worm.db?alt=media&token=96f46b2c-f21f-4d23-b51f-ece048ff896d)

## Preparation
If you need to make original dataset for other species or other data resource, prepare

1. Download rhe reference sequence  (.fa, .fasta). 


* Reference annotation
TAB (0x09) delimited text file

||||
|--|--|--|
||||



### Gene annotation


Prebuild plugin to load GFF format file and other data

Please refer the [plugin source]() for worm (C. elegans) as an example.



```sh
cd Moirei/src/plugin



```


