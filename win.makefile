
clean:
	-del /F /Q target

all: clean
	xcopy additions target
	$(MAKE) -C src
