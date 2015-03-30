plugins/color-chrs.so: plugins/color-chrs.c version.h version.c HMM.h HMM.c $(htslib_shared)
	$(CC) -fPIC -shared $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) -L$(HTSDIR) $(LDFLAGS) -o $@ HMM.c version.c $< -lhts $(LIBS)
