function []=ConsGaussFit05082016(filename,outfolder,outpdf,outcsv)
	%modify the removal of peaks with heights less than 1/5 of the max
%This function will create a folder that contains the plot of the optimal Gaussian fitted values and peakfinder. It will also save a table that contains the location and the height of the fitted peaks.
%added: plot component peaks.


% The format of the input data table:
% should have one string column of protein ID and a header row.

% 	CURVE FITTING ALGORITHM
% 1)	If a profile has less than 2 adjacent non-zero fractions, skip this profile.
% 2)	Decide the maximum number of Gaussian peaks: 3~5 non-zero fractions for up to 1 Gaussian peak, 6~8 for up to 2 peaks, 9~11 for 3 peaks and 12 or above for 4 peaks.
% 3)	Find the starting points:
% 	a)	Find fractions that have intensities not less than that of the left and right fraction and set them as "peak candidates"
% 	b)	Select the global maximum as the 1st starting point
% 	c)	Exclude fractions that are 3 or less fractions to the selected fractions and choose maximum from the rest of the "peak candidates".
% 	d)	Repeat step c until the number of peaks reaches the number decided in step 2.
% 4)	Fit several Gaussian curves, the maximum number of peaks is decided in step 2.
% 5)	Use the Gaussian curves that have lowest BIC (a combination of sum of squared error and the number of parameters) as the final choice of number of Gaussian curves.
% 6)	Remove the peaks that have intensities less than 1/5 of the maximum and show the final result.


	ProtTab=importdata(filename);
	allChromats=ProtTab.data;
	numChromats=size(ProtTab.data);numChromats=numChromats(1);
	ProtID=ProtTab.textdata(2:(numChromats+1),1);
	mkdir(outfolder);
	set(gcf, 'Visible', 'off')
	fileName1=strcat('./',outfolder,'/',outcsv,'.csv');
	fileName1=fopen(fileName1,'w');
	format='"%s", %.4g,%.4g,%.4g,%.4g, %g,%g, %g,%g,%g,%g, %.3g,%.3g,%.3g,%.3g, %.5g,%.5g,%.5g,%.5g, %.5g,%.5g,%.5g,%.5g, %.5g,%.5g,%.5g,%.5g\n';
	locHi={'ProtID' 'BIC1' 'BIC2' 'BIC3' 'BIC4' 'NumPeakLowBIC' 'NumPeakDelLowIntensity' 'Start1' 'Start2' 'Start3' 'Start4' 'Frac1' 'Frac2' 'Frac3' 'Frac4'  'Height1' 'Height2' 'Height3' 'Height4' 'Width1' 'Width2' 'Width3' 'Width4' 'Area1' 'Area2' 'Area3' 'Area4'};
	% Columns in the .csv table:
	% 1.	ProtID: protein ID
	% 2.	BIC1-4: the BIC of the fits of 1 to 4 fitted Gaussian curves
	% 3.	NumPeakLowBIC: the number of peaks that have lowest BIC
	% 4.	NumPeakDelLowIntensity: the number of peaks after deleting those which have intensity lower than the 1/5 of the maximum intensity
	% 5.	Start1-4: the starting peak location for Gaussian fitting
	% 6.	Frac1-4: the location of fitted Gaussian curves
	% 7.	Height1-4: the heights of the fitted Gaussian curves.
	% 8.	Width1-4: the widths of the fitted Gaussian curves.
	% 9.	Area1-4: the areas under the curves of the fitted Gaussian curves.

	fprintf(fileName1,'"%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s","%s"\r\n',locHi{:});
	%numChromats=10
	for chromatograph = 1:numChromats 
		y = allChromats(chromatograph,:);
		y=y';
		x = (1:numel(y))';
		Highcon = max( diff([0 (find( ~ (y'>0))) numel(y')+1])-1); %max number of consecutive positive
		if(Highcon<2)%must have at least 2 consecutive non-zero fractions to fit 1 Gaussian
			continue
		end
		%the bic combine the sum of squared errors and the penalty of the number of parameters
		bic1=0;
		bic2=0;
		bic3=0;
		bic4=0;
		n=numel(y);
		mH = max(y); %finding max height in chromatogram;
		mu=y'/sum(y)*x; %the estimate of the mean, the location of gaussian
		sd=sqrt(2)*sqrt((x-y'/sum(y)*x)'.^2*y/sum(y));%the estimate of the peak width
		yleft=[0 ;y(1:(n-1))];
		yright=[y(2:n) ;0];
		peakset=y>=yleft &y>=yright; %choose those fractions that are large than the left and right neighbor
		b1=find(y==max(y));b1=b1(1);
		delpeak=y;
		peakset(intersect(1:n,(b1-3):(b1+3)))=0;
		delpeak(~peakset)=0;
		peakall=peakfinder(y);
		rsq=repmat(0,1,numel(peakall));
		for peakloc = peakall'
			[f1,g1] = fit(x, y,(fittype('gauss1')),(fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0 1 0.71],'Upper',[mH n 5],'StartPoint',[mH peakloc sd]))); 
			rsq(find(peakloc == peakall))=g1.rsquare;
		end
		maxindex=find(rsq == max(rsq));
		[f1,g1] = fit(x, y,(fittype('gauss1')),(fitoptions('method','NonlinearLeastSquares','Robust','off','Lower',[0 1 0.71],'Upper',[mH n 5],'StartPoint',[mH peakall(maxindex(1) ) sd] )) );
		bic1=n*(log(2*pi)+log(g1.sse)-log(n))+n+3*log(n);
		sortf1=f1.b1;
		Npeak=1;
		if(max(delpeak)~=0 && sum(y>0)>=6) 
			%fit at least 2 Gaussians if number of non-zero fractions is at least 6
			b2=find(delpeak==max(delpeak));b2=b2(1);
			peakset(intersect(1:n,(b2-3):(b2+3)))=0;
			delpeak(~peakset)=0;  
			sorted=sort([b1 b2]);
			sb1=sorted(1);	sb2=sorted(2);
			sd1=max([sd/2 0.71]);sd2=max([sd/2 0.71]);
			ft=fittype('consExp(a1, a2, b1, b2, c1, c2, x)');
			[f2,g2]=fit(x, y,ft ,'Algorithm','Trust-Region', 'StartPoint',[y(sb1),y(sb2),sb1,sb2-4,sd1,sd2],...
				'Lower',[0,0, 1,1, 0.71,0.71],'Upper',[1.5*mH,1.5*mH,n-4,n-4,5,5]);
			bic2=n*(log(2*pi)+log(g2.sse)-log(n))+n+6*log(n);
			sortf2=sort([f2.b1,f2.b2])+[0 4];
			Npeak=2;
			if(max(delpeak)~=0 && sum(y>0)>=9)
				%fit at least 3 Gaussians if number of non-zero fractions is at least 9
				b3=find(delpeak==max(delpeak));b3=b3(1);
				peakset(intersect(1:n,(b2-3):(b2+3)))=0;
				delpeak(~peakset)=0;
				sorted=sort([b1 b2 b3]);
				sb1=sorted(1);	sb2=sorted(2);	sb3=sorted(3);
				sd1=max([sd/3 0.71]);sd2=max([sd/3 0.71]);sd3=max([sd/3 0.71]);
				ft=fittype('consExp3(a1, a2, a3, b1, b2, b3, c1, c2, c3, x)');
				[f3,g3]=fit(x, y,ft, 'Algorithm','Trust-Region', 'StartPoint',[y(sb1),y(sb2),y(sb3),sb1,sb2-4,sb3-8,sd1,sd2,sd3],...
					'Lower',[0,0,0, 1,1,1, 0.71,0.71,0.71],'Upper',[1.5*mH,1.5*mH,1.5*mH,n-8,n-8,n-8,5,5,5]);
				sortf3=sort([f3.b1,f3.b2,f3.b3])+[0 4 8];
				delpeak(intersect(1:n,(b3-3):(b3+3)))=0;
				bic3=n*(log(2*pi)+log(g3.sse)-log(n))+n+9*log(n);
				Npeak=3;
				if(max(delpeak)~=0 && sum(y>0)>=12)
					%fit at least 4 Gaussians if number of non-zero fractions is at least 12
					b4=find(delpeak==max(delpeak));b4=b4(1);
					sorted=sort([b1 b2 b3 b4]);
					sb1=sorted(1);	sb2=sorted(2);	sb3=sorted(3); sb4=sorted(4);
					sd1=max([sd/4 0.71]);sd2=max([sd/4 0.71]);sd3=max([sd/4 0.71]);sd4=max([sd/4 0.71]);
					ft=fittype('consExp4(a1, a2, a3, a4, b1, b2, b3, b4, c1, c2, c3, c4, x)');
					[f4,g4]=fit(x, y,ft, 'Algorithm','Trust-Region', 'StartPoint',[y(sb1),y(sb2),y(sb3),y(sb4),sb1,sb2-4,sb3-8,sb4-12,sd1,sd2,sd3,sd4],...
						'Lower',[0,0,0,0, 1,1,1,1, 0.71,0.71,0.71,0.71],'Upper',[1.5*mH,1.5*mH,1.5*mH,1.5*mH,n-12,n-12,n-12,n-12,5,5,5,5]);
					sortf4=sort([f4.b1,f4.b2,f4.b3,f4.b4])+[0 4 8 12];
					bic4=n*(log(2*pi)+log(g4.sse)-log(n))+n+12*log(n);
					Npeak=4;
				end
			end
		end
		bic=[bic1 bic2 bic3 bic4];
		bic=bic(bic>0);
		Npeak=find(bic==min(bic));
		plot(y,'k');	ylim([0,1.5*mH]);	xlim([0,n+2]);
		hold on
		switch Npeak
			case 1
				plot(f1,'g');
				legend(['raw data, starting:',num2str([b1])],['fitted Gauss locations:',num2str(sortf1,2)]);
				biccell=num2cell(bic);
				bloc=f1.b1;
				Hei=f1.a1;
				Wid=f1.c1;
				startloc=b1;
			case 2
				plot(f2,'g');
				legend(['raw data, starting:',num2str(sort([b1 b2]))],['fitted Gauss locations:',num2str(sortf2,2)]);
				[bloc ord]=sort([f2.b1 f2.b2]);
				bloc=bloc+[0 4];
				Hei=[f2.a1 f2.a2];
				Wid=[f2.c1 f2.c2];
				startloc=[b1 b2];
			case 3
				plot(f3,'g');
				legend(['raw data, starting:',num2str(sort([b1 b2 b3]))],['fitted Gauss locations:',num2str(sortf3,2)]);
				[bloc ord]=sort([f3.b1 f3.b2 f3.b3]);
				bloc=bloc+[0 4 8];
				Hei=[f3.a1 f3.a2 f3.a3];
				Wid=[f3.c1 f3.c2 f3.c3];
				startloc=[b1 b2 b3];
			
			case 4
				plot(f4,'g');
				legend(['raw data, starting:',num2str(sort([b1 b2 b3 b4]))],['fitted Gauss locations:',num2str(sortf4,2)]);
				[bloc ord]=sort([f4.b1 f4.b2 f4.b3 f4.b4]);
				bloc=bloc+[0 4 8 12];
				Hei=[f4.a1 f4.a2 f4.a3 f4.a4];
				Wid=[f4.c1 f4.c2 f4.c3 f4.c4];
				startloc=[b1 b2 b3 b4];
		end	
	
		bloc=bloc(Hei>mH/5);
		Wid=Wid(Hei>mH/5);
		Hei=Hei(Hei>mH/5);
		%remove peaks lower than 1/5 of maximum intensity
		Area=normprob(n+0.5,0.5,Hei,bloc,Wid);
		Area=num2cell(Area);
		Hei=num2cell(Hei);
		Wid=num2cell(Wid);
		biccell=num2cell(bic);			
		locHi=[ProtID(chromatograph) biccell cell(1,4-numel(bic)) [Npeak] [numel(Hei)] num2cell(sort(startloc)) cell(1,4-numel(startloc)) num2cell(bloc) cell(1,4-numel(Hei)) Hei cell(1,4-numel(Hei))  Wid cell(1,4-numel(Wid)) Area cell(1,4-numel(Area))];
		fprintf(fileName1,format,locHi{:});

		title([strcat('Optimal ',num2str(Npeak),' Fitted Gauss, protein',ProtID(chromatograph)),strcat('bic: ',num2str(bic,3)),strcat('peaks after deleting less than one fifth of global max: ',num2str(bloc,2))]);
		hold off
		set(gca, 'LooseInset', get(gca, 'TightInset'));
		set(gcf, 'PaperUnits', 'inches');
		saveas(gcf, [strcat('./',outfolder,'/'),num2str(chromatograph),'plotAllGaussOpt'], 'pdf');

		Hei=cell2mat(Hei);
		Wid=cell2mat(Wid);
		x=(0:n*10)/10;
		for i = 1: numel(Hei)
			plot(y,'k');	ylim([0,1.5*mH]);	xlim([0,n+2]); hold on;
			plot(x, Hei(i)*exp(-((x-bloc(i))/(Wid(i))).^2),'g');
			Area=normprob(n+0.5,0.5,Hei(i),bloc(i),Wid(i));
			title([strcat('Optimal number of peaks: ',num2str(numel(Hei)),' Component Gaussian peak of protein',ProtID(chromatograph)), strcat('Component number:',num2str(i)),strcat('Area Under the curve',num2str(Area,4))]);
			hold off
			set(gca, 'LooseInset', get(gca, 'TightInset'));
			set(gcf, 'PaperUnits', 'inches');
			saveas(gcf, [strcat('./',outfolder,'/'),num2str(chromatograph),'plotAllGaussOpt',num2str(i)], 'pdf');
		end
	end

	j=1;
	f=cell(1,1);
	for i=1:numChromats
		if(exist(strcat('./',outfolder,'/',num2str(i),'plotAllGaussOpt.pdf'))==2)
			f{j}=strcat('./',outfolder,'/',num2str(i),'plotAllGaussOpt.pdf');
			j=j+1;
		end
	end
	fclose('all');
	append_pdfs(strcat('./',outfolder,'/',outpdf,'.pdf'), f{:});
	delete(strcat('./',outfolder,'/*plotAllGaussOpt.pdf'))

	j=1;
	f=cell(1,1);
	for i=1:numChromats
		if(exist(strcat('./',outfolder,'/',num2str(i),'plotAllGaussOpt2.pdf'))~=2)
			continue
		end
		for k=1:4
			if(exist(strcat('./',outfolder,'/',num2str(i),'plotAllGaussOpt',num2str(k),'.pdf'))==2)
			f{j}=strcat('./',outfolder,'/',num2str(i),'plotAllGaussOpt',num2str(k),'.pdf');
			j=j+1;
			end
		end
	end
	fclose('all');
	append_pdfs(strcat('./',outfolder,'/',outpdf,'ComponentPeak.pdf'), f{:});
	delete(strcat('./',outfolder,'/*plotAllGaussOpt*.pdf'))



end
