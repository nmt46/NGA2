function st = read_stats(filepath)
% Read the contents of a statistics save file and return the reconstructed
% stats object
% Inputs:
%   filepath = character array of file path to the statistics file
% Outputs:
%   st = stats parent object reconstructed from statistics file

myF = fopen(filepath);

% Read stats header
st.name = strtrim(convertCharsToStrings(char(fread(myF,64,'char').')));
st.nLoc =     fread(myF,1,'int32');
st.nAr  =     fread(myF,1,'int32');
st.nDef =     fread(myF,1,'int32');
st.sum_time = fread(myF,1,'double');

% Loop over definitions
for iDef = 1:st.nDef
    st.def(iDef).name = strtrim(convertCharsToStrings(char(fread(myF,64,'char').')));
    st.def(iDef).expP = fread(myF,6,'int32');
    st.def(iDef).nTA  = fread(myF,1,'int32');
    st.def(iDef).iAs  = fread(myF,st.def(iDef).nTA,'int32');
end

% Loop over array pointers and get names
for n = 1:st.nAr
    st.arp(n).name = strtrim(convertCharsToStrings(char(fread(myF,64,'char').')));
end

% Loop over stations
for iLoc = 1:st.nLoc
    st.stats(iLoc).name = strtrim(convertCharsToStrings(char(fread(myF,64,'char').')));
    st.stats(iLoc).dim  = fread(myF,1,'int32');
    st.stats(iLoc).imin = fread(myF,1,'int32');
    st.stats(iLoc).imax = fread(myF,1,'int32');
    st.stats(iLoc).jmin = fread(myF,1,'int32');
    st.stats(iLoc).jmax = fread(myF,1,'int32');
    st.stats(iLoc).kmin = fread(myF,1,'int32');
    st.stats(iLoc).kmax = fread(myF,1,'int32');
    dims = fread(myF,3,'int32');
    st.stats(iLoc).vals = zeros(dims.');
    for n = 1:st.nDef
        for k = 1:dims(3)
            for j = 1:dims(2)
                st.stats(iLoc).vals(:,j,k,n) = fread(myF,dims(1),'double');
            end
        end
    end
end

fclose(myF);
end