Like saanity check but for mean field






function sanity2 = sanitychkmf(jn,Sstruct,JHstruct)
   sanity2 = struct('th',{},'tchk',{},'mtchk',{},'mimj',{},'chi',{},'mchi',{},'saneh',{},'sanechi',{},'Jmk',{},'mi',{});
   for i = 1:jn
    	h = JHstruct(i).Hsparse;
        J = JHstruct(i).Jsparse;
    	mk = Sstruct(i).mfinal;
        
    	sisj = Sstruct(i).Cfinal;
    	
        tchk = tanh(h) - mk; 
        mtchk = mean(tchk);
        
        Jmk = mk*J;
        mi = tanh(h+Jmk);

        mimj = mi'*mi;
    	chi = sisj - mimj;
        
    	mchi = mean(chi,'all');

        if abs(mtchk) < 0.01
            saneh = 1;
        else
            saneh = 0;
        end

    	if abs(mchi) < 0.01
    		sanechi = 1;
    	else
    		sanechi = 0;
    	end

    	sanity2(i).th = tanh(h);
        sanity2(i).tchk = tchk;
        sanity2(i).mtchk = mtchk;
    	sanity2(i).mimj = mimj;
        sanity2(i).mi = mi;
        sanity2(i).Jmk = Jmk;
    	sanity2(i).chi = chi;
    	sanity2(i).mchi  = mchi;
        sanity2(i).saneh = saneh;
    	sanity2(i).sanechi = sanechi;

    end
