function data_dsorted=sort_data_by_depth(data,depth_of_data, ...
                                          deps_of_interest, binwidth)
 % Description:
 % data_dsorted=sort_data_by_depth(data,depth_of_data,deps_of_interest,
 % binwidth);
 % Input: 
 %   1. data: data to be sorted by depth (a vector)
 %   2. depth_of_data: the associate depth for each data entry.
 %   3. deps_of_interest: the depths at which the data is sorted by
 %   4. binwidth:  depth between deps_of_interest +- binwidth is selected.
 % Output:
 % data_dsorted: store in a cell array.
                                      
for i=1:length(deps_of_interest)
    current_dep=deps_of_interest(i);
    dep_dif=abs(depth_of_data - current_dep);
    depID=find(dep_dif <= round(deps_of_interest(i)*binwidth));
    dep_top=deps_of_interest(i)*(1+binwidth);
    dep_bot=deps_of_interest(i)*(1-binwidth);
    if ~isempty(depID)
        prctg=length(depID)/length(depth_of_data) * 100;
        disp(['found' num2str(prctg,'%4.1f')  'of data between' ...
            num2str(dep_bot,'%5.1f') '~' num2str(dep_top,'%5.1f')]);
    else
        disp(['no data is found between'  num2str(dep_bot,'%5.1f') '~' num2str(dep_top,'%5.1f')]);
    end
    
    % use this ID vector to select data out of each variable:
    data_dsorted{i}=data(depID);
    
end

end