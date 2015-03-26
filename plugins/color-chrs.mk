plugins/color-chrs.so: plugins/color-chrs.c version.h version.c HMM.h HMM.c $(htslib_shared)
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ HMM.c version.c $< -L$(HTSDIR) -lhts
