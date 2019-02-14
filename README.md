
[<img align="center" src="NBIS.png" width="200" height="100" />](https://nbis.se) 
<h2><em>G</em>enome <em>A</em>ssembly <em>A</em>nnotation <i>S</i>ervice (GAAS)</h2>  
Contains development done in the GAA (Genome Assembly Annotation) Service.

---------------------------

## [__annotation__](annotation)  
Annotation directory contains development related to annotation side of the service.  

## [__assembly__](assembly)  
Assembly directory contains development related to assembly side of the service.  

---------------------------

#### Installation

  * Clone and install GAAS
  
To use the tools in this repository, clone the directory and run `make install` to update some paths:
```
git clone https://github.com/NBISweden/GAAS.git
cd GAAS
```

  * A) For the use through ***env***  

    * A.1) Updates paths in the environment profiles to point to the correct GAAS repository location  
    ```
    make install   
    ```

    * A.2) Dependencies
      *  You might check that all dependencies are filled up. Depending the scripts you want to use, all dependencies are not required.  
      ```
      make check
      ```
      * install the missing dependencies of your choice

    * A.3) Load the correct profiles (add NBIS libs and tools to the PATH)  
    Three profiles are available to setup the necessary environment variables to use the scripts:

      * By default:
      ```
      source profiles/activate_env
      ```

       * If you are on Rackham:
      ```
      source profiles/activate_rackham_env
      ```

       * If you are on NBIS's servers:
      ```
      source profiles/activate_nbis_env
      ```
    * A.3) To get out of the nbis environment and restore your previous environment type  
  
     ```
     deativate
     ```

  * B) For a permanent use  
  
      * B.1) You might check that all dependencies are filled up. Depending the scripts you want to use, all dependencies are not required.
      ```
      make check
      ```
      
      * B.2) Add the path to the BILS perl library as well as the bin folder containing all tools. You can add in you *~/.bashrc* or *~/.profile* file.
      ```
      export PERL5LIB=$PERL5LIB:/pathTo/GAAS/annotation
      export PATH=${PATH}:/pathTo/GAAS/annotation/Tools/bin
      ```
      
      
