function plotruntimevsdim(dmin, dmax, time)
    mtime=mean(time,1);
    plot((dmin:dmax),mtime,'black')
    hold on
    ylabel('runtime (seconds)')
    xlabel('dimension')
end