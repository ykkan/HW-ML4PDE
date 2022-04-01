function ploterrorvsruntime(v,value,time)
    merror=mean(abs(value-v))/abs(v);
    mtime=mean(time);
    loglog(mtime,merror,'black-o');
    hold on
    loglog(mtime,1./(mtime).^(1/3)*mtime(1)^(1/3)*merror(1),'black');
    ylabel('relative approximation error');
    xlabel('runtime (seconds)');
    legend('relative approximation error','slope -1/3');
end
