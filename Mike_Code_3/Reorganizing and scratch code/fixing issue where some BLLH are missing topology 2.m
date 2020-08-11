%% For fixing issue where some BLLH are missing topology 2 ()
for nnn = 61: length(multiallstruct);
	if length(multiallstruct(nnn).BLLH) == 5
		multiallstruct(nnn).BLLH = [multiallstruct(nnn).BLLH(1),multiallstruct(nnn).BLLH(1),multiallstruct(nnn).BLLH([2:end])]
	end
end
