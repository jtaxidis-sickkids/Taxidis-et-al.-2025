function [animal,folder,day] = get_animal_name(day)

[parent_folder,day] = fileparts(day);
folder = fileparts(parent_folder);
[~,animal] = fileparts(folder);

k = strfind(day,'_');       % There should be 2 underscores
if length(k) > 2
    day(k(3):end) = [];
end