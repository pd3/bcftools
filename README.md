This branch contains experimental version of bcftools mpileup and gVCF calling.

```
git clone --branch=exp/gvcf --recursive git://github.com/pd3/bcftools.git
cd bcftools; make 

cd bcftools
git pull && git submodule update --recursive
make clean && make 
```

