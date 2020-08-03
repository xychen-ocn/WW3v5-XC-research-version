function [outlierVal,idx]=locate_outliers_in_data(data,criterion1, criterion2)

switch nargin
    case 2
        % find data larger than criteria:
        idx=find(data>1.05*criterion1);
        if ~isempty(idx)
            outlierVal=data(idx);
        else
            outlierVal=[];
            idx=[];
        end
    case 3
        % find data larger than criteria:
        id1=find(data>1.05*criterion1);
        if ~isempty(id1)
            outlierVal1=data(id1);
        else
            outlierVal1=[];
            id1=[];
        end
        
        % find data smaller than criteria:
        id2=find(data<0.9*criterion2);
        if ~isempty(id2)
            outlierVal2=data(id2);
        else
            outlierVal2=[];
            id2=[];
        end
        
        idx=[id1 id2];
        outlierVal=[outlierVal1 outlierVal2];
        
    otherwise
        disp('not enough input arguements!');
end
    

return