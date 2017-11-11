function MiddleAxis(ticksize,tick,range,offsetx,offsety)
    line(range(1:2),[0 0],'color','k')
    for i = 2 :length(tick)
        line([tick(i) tick(i)],[-ticksize ticksize],'color','k');
        text(tick(i)-offsetx,offsety,num2str(tick(i)));
    end

end