function data=getPatientDataFromStruct(data,fieldname,searchkey)
    %getPatientData - searches a struct array for a particular patient data
    %
    %PARAMS:
    %data - the struct array to search
    %fieldname - the field variable to search
    %searchkey - the value to search for
    %
    %EXAMPLE:
    %patients(1).name='John Doe';
    %patients(1).age=52;
    %patients(2).name='Jane Doe';
    %patients(2).age=12;
    %
    %ptn=getPatientData(patients,'name','Jane Doe')
    %ptn.age
    
    searchVals = {data(:).(filedname)};
    idx = find(ismember(searchVals,'searchkey'));
    data = data(idx);