# svn $Id$
#::::::::::::::::::::::::::::::::::::::::::::::::::::: Hernan G. Arango :::
#   Copyright (c) 2002-2013 The ROMS/TOMS Group                         :::
#   Licensed under a MIT/X style license                                :::
#   See License_ROMS.txt                                                :::
#::::::::::::::::::::::::::::::::::::::::::::: Benjamin K. Galton-Fenzi :::

local_sub  := ROMS/Nonlinear/IceShelf

local_lib  := libNLM_is.a
local_src  := $(wildcard $(local_sub)/*.F)

$(eval $(call make-library,$(local_lib),$(local_src)))

$(eval $(compile-rules))
