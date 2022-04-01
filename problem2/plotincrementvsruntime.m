function plotincrementvsruntime(value,time)
    [~,rhomax]=size(value);
    merror=mean(abs(value(:,2:rhomax)-value(:,1:rhomax-1)))/abs(mean(value(:,rhomax)));
    mtime=mean(time(:,1:rhomax-1));
    loglog(mtime,merror,'black-o');
    hold on
    loglog(mtime,1./(mtime).^(1/3)*mtime(1)^(1/3)*merror(1),'black');
    ylabel('relative approximation increments');
    xlabel('runtime (seconds)');
    legend('relative approximation increments','slope -1/3');
end
