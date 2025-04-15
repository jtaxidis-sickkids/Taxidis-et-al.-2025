function [sessiontypes, Nprobes] = Set_Sessiontypes(session)

%% SESSIONTYPES
switch session    
    % Naive No-Opto
    case {'ZD028_0822','ZD029_0823','ZD028_0824','ZD029_0824'}
        sessiontypes = repmat("no light",[8,1]);
    % -----------------------------
        
    % Naive Opto
    case {'CD235_0825','CD234_0826','CD235_0826','ZD031_0827','ZD032_0827','ZD031_0828','ZD032_0828'}
        sessiontypes = {'no light','no light','no light','no light','no light',...
            'rebound','rebound','rebound','rebound',...
            'onset','onset','onset','onset',...
            'full odor','full odor','full odor'}';
        
    case 'CD234_0825'
        sessiontypes = {'no light','no light','no light','no light','no light',...
            'full odor','full odor','full odor',...
            'rebound','rebound','rebound','rebound',...
            'onset','onset','onset','onset'}';
    % -----------------------------
        
    % Trained No-Opto
    case 'ZD028_0904'
        sessiontypes = repmat("no light",[11,1]);
        
    case 'ZD029_0904'
        sessiontypes = repmat("no light",[16,1]);
    % -----------------------------
        
    % Trained Opto
    case {'CD234_0905','CD235_0905','ZD031_0906'}
        sessiontypes = {'no light',...
            'rebound','rebound', 'onset','onset', 'full odor','full odor',...
            'rebound','onset','full odor',...
            'no light'}';
        
    case {'CD235_0906','ZD031_0907','CD234_0907','ZD031_0908'}
        sessiontypes = {'no light',...
            'rebound','rebound', 'onset','onset', 'full odor','full odor',...
            'hyper','hyper',...
            'rebound','onset'}';
    % -----------------------------
end

%% NPROBES
switch session    
    case {'ZD028_0904','ZD029_0904','CD234_0905','ZD031_0906','CD235_0906'}
        Nprobes = 2;
    otherwise
        Nprobes = 1;
end
