
#
#   make - Names
#

    MAKE_NAME=libinter

#
#   make - Directories
#

    MAKE_BINARY=bin
    MAKE_DOCUME=doc
    MAKE_OBJECT=obj
    MAKE_SOURCE=src

#
#   make - Files
#

    MAKE_SCC=$(wildcard $(MAKE_SOURCE)/*.c)
    MAKE_OBC=$(addprefix $(MAKE_OBJECT)/, $(addsuffix .o, $(notdir $(basename $(MAKE_SCC)))))

#
#   make - Constants
#

    MAKE_CCC=gcc
    MAKE_LKD=ar
    MAKE_OPC=-Wall -c -funsigned-char -std=c99
    MAKE_DOC=doxygen

#
#   make - All
#

    all:directories $(MAKE_NAME)

#
#   make - Binaries
#

    $(MAKE_NAME):$(MAKE_OBC) $(MAKE_OBP)
	$(MAKE_LKD) rcs $(MAKE_BINARY)/$(MAKE_NAME).a $^

#
#   make - Objects
#

    $(MAKE_OBJECT)/%.o:$(MAKE_SOURCE)/%.c
	$(MAKE_CCC) $(MAKE_OPC) -o $@ $<

#
#   make - Documentation
#

    documentation:directories clean-doc
	rm $(MAKE_DOCUME)/html -rf && $(MAKE_DOC)

#
#   make - Directories
#

    directories:
	mkdir -p $(MAKE_BINARY)
	mkdir -p $(MAKE_DOCUME)
	mkdir -p $(MAKE_OBJECT)

#
#   make - Clean
#

    clean:
	rm $(MAKE_BINARY)/* -f
	rm $(MAKE_OBJECT)/*.o -f

