# make file for Habitats

CC = g++
WARN = -Wall
OPTI = -O0
DBTOK = -g
CFLAGS = -I/usr/local/lib/R/include -I/usr/local/lib64/R/include -fPIC -ffast-math 
pkgName = ecomix

ippmDIR = /home/woo457/Dropbox/ecomix/developing_sam_src

$(pkgName).so: sam_ippm.o sam_ippm_data.o sam_ippm_fits.o sam_ippm_params.o sam_ippm_derivs.o sam_ippm_all_classes.o sam_ippm_control.o  
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -shared -o $(pkgName).so sam_ippm.o sam_ippm_data.o sam_ippm_fits.o sam_ippm_params.o sam_ippm_derivs.o sam_ippm_all_classes.o sam_ippm_control.o 

sam_ippm_data.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm_data.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c $(ippmDIR)/sam_ippm_data.cc

sam_ippm_control.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm_control.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c $(ippmDIR)/sam_ippm_control.cc

sam_ippm_fits.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm_fits.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c $(ippmDIR)/sam_ippm_fits.cc

sam_ippm_params.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm_params.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c  $(ippmDIR)/sam_ippm_params.cc

sam_ippm_derivs.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm_derivs.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c  $(ippmDIR)/sam_ippm_derivs.cc

sam_ippm_all_classes.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm_all_classes.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c  $(ippmDIR)/sam_ippm_all_classes.cc

sam_ippm.o: $(ippmDIR)/sam_ippm.h $(ippmDIR)/sam_ippm.cc
	$(CC) $(CFLAGS) $(OPTI) $(DBTOK) $(WARN) -c $(ippmDIR)/sam_ippm.cc

clean:
	rm -f  $(ippmDIR)/*.o $(ippmDIR)/*.so *.o *.so
