

latency=[EEG.event.latency];
pre=latency(find(ismember({EEG.event.type},{'4' '5'})==1));
pst=latency(find(ismember({EEG.event.type},{'4' '5'})==1)+1);
rt=pst-pre;