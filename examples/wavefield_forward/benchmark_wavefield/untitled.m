% Opens editor's state file
    filepart = sprintf('MathWorks\\MATLAB\\R%s\\%s', version('-release'), 'MATLAB_Editor_State.xml');
    filename = fullfile(getenv('APPDATA'), filepart);
    document = xmlread(filename);

    % Get information about 'File' nodes
    recentFiles = struct([]);
    fileNodes = document.getElementsByTagName('File');
    for fni = 1:(fileNodes.getLength())

       attributes = fileNodes.item(fni-1).getAttributes(); % Careful, zero based indexing !

       for ai = 1:(attributes.getLength())

           % Get node attribute
           name = char(attributes.item(ai-1).getName()); % Zero based + need marshaling COM 'string' type
           value = char(attributes.item(ai-1).getValue()); % Zero based + need marshaling COM 'string' type

           % Save in structure
           name(1) = upper(name(1)); % Just because I prefer capital letter for field names ...
           recentFiles(fni).(name) = value;

       end

    end