%% to add symmetery averaging and rerr to multiallstruct for putting in figure
%%



function multiallstruct = addsymm(multiallstruct)

models = [{'BIC'     },{'AIC'     },{'MDLl'    },{'MDLu'    },{'MDLent'  },{'MDLcount'}];

topologies = [1,3,4,5,6];


%for i = 1:length(multiallstruct)
for i = 1:length(multiallstruct) % special case where most of the struct is already done

	%for t =1:length(topologies)	%topology should only be 5 since we are skipping the topo=2 atm
	for t =1:length(topologies)		% one above fails to when querying for 
	Jtru = multiallstruct(i).Jcontru{t};
	topo = topologies(t);
	
		for m =1:length(models)	

			modelname = models{m};

			%%need to fix the asymmetry finder
			%------------------------------------------------------------
			row = multiallstruct(i).BLLH(topo).asymrow.(modelname) ;
			col = multiallstruct(i).BLLH(topo).asymcol.(modelname) ;

			Jrecon 	=	multiallstruct(i).BLLH(topo).Jrecon.(modelname)  ;
			Jsymrecon = zeros(size(Jrecon));
			N = size(Jrecon,1);
			
			for coordno = 1:length(row)
            	ci = row(coordno);
            	
				wml = multiallstruct(i).BLLH(topo).NodeModel(ci).w_ML(1:N,2:end) ;
				
				idx0 = 1:1:size(wml,2);
				zeron = zeros(size(wml,1),1);
				A = setdiff(idx0,ci:N);
				B = setdiff(idx0,1:ci-1);

				wml2 = [wml(:,A),zeron,wml(:,B)];
				
	
				v = nonzeros(wml2(:,col(coordno)));
				reconval = v(end);
	
				Jsymrecon(row(coordno),col(coordno)) = reconval;
			end

			multiallstruct(i).BLLH(topo).Jsymrecon.(modelname) = Jsymrecon + Jrecon;
			 Jsymcon = (Jsymrecon + Jrecon) ~= 0 ;
			 multiallstruct(i).BLLH(topo).Jsymcon.(modelname) =Jsymcon


			%------------------------------------------------------------
				
			%Jtru = multiallstruct(i).Jcontru{t};
			
%			Jrecon 	=	multiallstruct(i).BLLH(topo).Jrecon.(modelname)  ; % moving above
			
			Jsymrecon =	multiallstruct(i).BLLH(topo).Jsymrecon.(modelname);
			
			%% For making the avg symmetry Jrecon
			avgcon = (triu(Jsymrecon) + tril(Jsymrecon)')/2;
			Javgsymrecon = avgcon + avgcon';
			symtotconerr = sum(double(not((Javgsymrecon - Jtru) == 0)),'all');
			symperconerr = symtotconerr ./ N.^2;
			
			multiallstruct(i).BLLH(topo).Javgsymrecon.(modelname) = Javgsymrecon;
			multiallstruct(i).BLLH(topo).symtotconerr.(modelname)      = symtotconerr;
			multiallstruct(i).BLLH(topo).symperconerr.(modelname)      = symperconerr;


			multiallstruct(i).BLLH(topo).Jrrerr.(modelname) 		=sqrt(sum((Jrecon 		- Jtru).^2 / sum(Jtru.^2)));
			multiallstruct(i).BLLH(topo).Jsymrrerr.(modelname) 	=sqrt(sum((Jsymrecon 		- Jtru).^2 / sum(Jtru.^2)));
			multiallstruct(i).BLLH(topo).Javgsymrrerr.(modelname) =sqrt(sum((Javgsymrecon 	- Jtru).^2 / sum(Jtru.^2)));
			
		
		
end	
end
end