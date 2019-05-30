function VerifyVersion(RequestedVersion)
    
    try
    Version = Tools.Version();
    if ~strcmp(Version,RequestedVersion)
        error('MDP:wrong version of Tools, expected version %s, found version %s',RequestedVersion,Version);
    end
catch err
    if strcmp(err.identifier, 'MATLAB:undefinedVarOrClass')
        error('MDP: please add parent folder to the path');
    else%if strcmp(err.identifier,'MDP:wrong version of Tools')
        sprintf(err.message);
        rethrow(err);
    end
end
    
end