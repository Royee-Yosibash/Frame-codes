close all;
fontSize = 26;

for iResults=1:numel(TotResults)
    h = figure();   haxes = axes(h); hold(haxes, 'on'); grid(haxes, 'on');
    for iCode = 1:numel(codeTypes)
        TotResults(iResults).plotByMeasure('MSE', codeTypes{iCode}, haxes);
    end
    legend(haxes);
    set(haxes, 'YScale', 'log');
    haxes.FontSize = fontSize;
end

for iResults=1:numel(TotResults)
    h = figure();   haxes = axes(h); hold(haxes, 'on'); grid(haxes, 'on');
    for iCode = 1:numel(codeTypes)
        TotResults(iResults).plotByMeasure('Frobenius norm', codeTypes{iCode}, haxes);
    end
    legend(haxes);
    set(haxes, 'YScale', 'log');
    haxes.FontSize = fontSize;
end

for iResults=1:numel(TotResults)
    h = figure();   haxes = axes(h); hold(haxes, 'on'); grid(haxes, 'on');
    for iCode = 1:numel(codeTypes)
        TotResults(iResults).plotByMeasure('mean(conditionNumber)', codeTypes{iCode}, haxes);
        TotResults(iResults).plotByMeasure('max(conditionNumber)', codeTypes{iCode}, haxes);
        TotResults(iResults).plotByMeasure('min(conditionNumber)', codeTypes{iCode}, haxes);
    end
    legend(haxes);
    set(haxes, 'YScale', 'log');
    haxes.FontSize = fontSize;
end
