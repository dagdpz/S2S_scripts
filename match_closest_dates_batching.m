function [out_files tolerated_indexes]=match_closest_dates_batching(exp_files,bas_files,tolerance)

for b=1:numel(bas_files)
rarranged_bas(b)=datenum(bas_files{b}{1}(end-7:end),'YYYYmmdd');
end

for e=1:numel(exp_files)
    date=datenum(exp_files{e}{1}(end-7:end),'YYYYmmdd');
    
    [cm,idx]=min(abs(rarranged_bas-date));
    out_files(e)=bas_files(idx);
    if cm>tolerance
      tolerated_indexes(e)=false;  
    else
      tolerated_indexes(e)=true; 
    end
end
end