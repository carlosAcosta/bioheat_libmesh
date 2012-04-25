
# $Id: Makefile,v 1.9 2008-10-16 18:04:01 roystgnr Exp $

# The location of the mesh library
LIBMESH_DIR ?= ../..
#LIBMESH_DIR ?= ../../../

# include the library options determined by configure
include $(LIBMESH_DIR)/Make.common

target 	   := ./bioheat-$(METHOD)

# include my common application options
include Make_template

