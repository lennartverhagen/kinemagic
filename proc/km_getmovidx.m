function varargout = km_getmovidx(cfg,tseries)
%--------------------------------------------------------------------------
%
%
% This file is part of the KineMagic toolbox
% Copyright (C) 2010, Lennart Verhagen
% L.Verhagen@donders.ru.nl
% version 2010-01-01
%--------------------------------------------------------------------------

movnr	= 1;
marker	= {};
axis	= {};

if nargin > 1
    if isfield(cfg,tseries)
        cfg = cfg.(tseries);
    % legacy reasons: app > apt
    elseif strcmpi(tseries,'gripapp') && isfield(cfg,'gripapt')
        cfg = cfg.gripapt;
    elseif strcmpi(tseries,'gripapt') && isfield(cfg,'gripapp')
        cfg = cfg.gripapp;
    elseif strcmpi(tseries,'gripappvel') && isfield(cfg,'gripaptvel')
        cfg = cfg.gripaptvel;
    elseif strcmpi(tseries,'gripaptvel') && isfield(cfg,'gripappvel')
        cfg = cfg.gripappvel;
    elseif strcmpi(tseries,'gripappacc') && isfield(cfg,'gripaptacc')
        cfg = cfg.gripaptacc;
    elseif strcmpi(tseries,'gripapacct') && isfield(cfg,'gripappacc')
        cfg = cfg.gripappacc;
    end
end
if isfield(cfg,'movnr'),	movnr	= cfg.movnr;	end
if isfield(cfg,'marker'),	marker	= cfg.marker;	end
if isfield(cfg,'axis'),     axis	= cfg.axis;     end

if ~iscell(marker), marker  = {marker}; end
if ~iscell(axis),   axis    = {axis};   end

varargout = {movnr,marker,axis};