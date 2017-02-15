% --new code--
% --Time-stretching for stereo sound--
% 8/11/2014

[a,fs,bits]=wavread('DaDeMo_-_Grand_Piano_Fazioli_Major_Chords_Middle_Pitch.wav');    %read wave, set some basic parametrics
N=2048;
fftsize=N;
TimeStretchRatio=1;
k=8;
new_hopsize=N/k;
hopsize=round(new_hopsize/TimeStretchRatio);

w=hann(N,'periodic');
COLA_ratio=sum(w.^2)/new_hopsize;       
w=w/sqrt(COLA_ratio); 

row=size(a,1);

%----------------------------left channel---------------------------------%

%preparation for sound resynthesis
newsound_pointnumber=floor(TimeStretchRatio*row);
newsound_l=zeros(newsound_pointnumber,1);

pre_phase_ana=zeros(fftsize/2,1);
pre_phase_syn=zeros(fftsize/2,1);
unwrapdata=2*pi*hopsize*(0:(fftsize/2-1))'/fftsize;

pin=0;
pout=0;

firsttime=true;

while(pin+N<row)&&(pout+fftsize<newsound_pointnumber)
    
    yl=a((pin+1:pin+N),1);   %get a frame
    
    %window the buffer
    aw=yl.*w;
    
    %fft-shift
    afs=fftshift(aw);
    
    %fft
    af=fft2(afs);
    af=af(1:fftsize/2);
    
    %calaulate magnitude and phase
    realpart=real(af);
    imagpart=imag(af);
    mag=sqrt((realpart.^2)+(imagpart.^2));                        
    phase_ana=atan2(imagpart,realpart);
          
%---------------------------processing-------------------------------------

    if firsttime
        phase_syn=TimeStretchRatio.*phase_ana;
        firsttime=false;
    else
        phase_increment=(phase_ana-pre_phase_ana)-unwrapdata;
        principal_determination=mod(phase_increment+pi,2*pi)-pi;
        partials_freq=principal_determination/hopsize+(unwrapdata/hopsize);
        
        % Update the phase in each bin
        phase_syn=pre_phase_syn+new_hopsize*partials_freq;
    end;

    % Compute DFT of the synthesis frame
    af= mag.* exp(1j*phase_syn);

    % Remember phases
    pre_phase_ana=phase_ana;
    pre_phase_syn=phase_syn;
        
    
%---------------------------resynthesis------------------------------------

    af(fftsize/2+2:fftsize)=fliplr(af(2:fftsize/2)');
    
    %ifft
    aif=real(ifft2(af));    %why we must use 'real'?

    %fftshift
    new_afs=fftshift(aif);

    %window the buffer
    new_aw=new_afs.*w;
    
    %compose the newsound    
    newsound_l(pout+1:pout+fftsize)=newsound_l(pout+1:pout+fftsize)+new_aw;
    pin=pin+hopsize;
    pout=pout+new_hopsize;
    
end;





%---------------------------right channel---------------------------------%

%preparation for sound resynthesis
newsound_r=zeros(newsound_pointnumber,1);

pre_phase_ana=zeros(fftsize/2,1);
pre_phase_syn=zeros(fftsize/2,1);
unwrapdata=2*pi*hopsize*(0:(fftsize/2-1))'/fftsize;

pin=0;
pout=0;

firsttime=true;

while(pin+N<row)&&(pout+fftsize<newsound_pointnumber)
    
    yr=a((pin+1:pin+N),2);   %get a frame
    
    %window the buffer
    aw=yr.*w;
    
    %fft-shift
    afs=fftshift(aw);
    
    %fft
    af=fft2(afs);
    af=af(1:fftsize/2);
    
    %calaulate magnitude and phase
    realpart=real(af);
    imagpart=imag(af);
    mag=sqrt((realpart.^2)+(imagpart.^2));                        
    phase_ana=atan2(imagpart,realpart);
    %phase_ana=angle(af);
          
%---------------------------processing-------------------------------------

    if firsttime
        phase_syn=TimeStretchRatio.*phase_ana;
        firsttime=false;
    else
        phase_increment=(phase_ana-pre_phase_ana)-unwrapdata;
        principal_determination=mod(phase_increment+pi,2*pi)-pi;
        partials_freq=principal_determination/hopsize+(unwrapdata/hopsize);
        
        % Update the phase in each bin
        phase_syn=pre_phase_syn+new_hopsize*partials_freq;
    end;
    
    % Compute DFT of the synthesis frame
    af= mag.* exp(1j*phase_syn);

    % Remember phases
    pre_phase_ana=phase_ana;
    pre_phase_syn=phase_syn;
        
    
%---------------------------resynthesis------------------------------------

    af(fftsize/2+2:fftsize)=fliplr(af(2:fftsize/2)');
    
    %ifft
    aif=real(ifft2(af));    %why we must use 'real'?

    %fftshift
    new_afs=fftshift(aif);

    %window the buffer
    new_aw=new_afs.*w;
    
    %compose the newsound    
    newsound_r(pout+1:pout+fftsize)=newsound_r(pout+1:pout+fftsize)+new_aw;
    pin=pin+hopsize;
    pout=pout+new_hopsize;
    
end;

newsound=[newsound_l,newsound_r];
wavwrite(newsound,fs,bits,'output.wav');