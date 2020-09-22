function [reg_label,date,Ibar,Rbar,totCases] = data_read_dpc_regioni(path_folder)

%
%   [reg_label,date,Ibar,Rbar] = data_read_dpc_regioni(path_folder)
%
%   data_read_dpc_regioni legge i dati della protezione civile per
%   le regioni
%
%   INPUTS:
%   path_folder : percorso folder dove ci sono i csv
%
%   OUTPUTS:
%   reg_label   : etichette regioni
%   date        : vettore stringhe di date
%   Ibar        : cella di dimensioni(21,1) dove in {i,:}(:) ho Ibar della
%                 regione i
%                 totale casi positivi, I(t) secondo il SIR
%   Rbar        : cella di dimensioni(21,1) dove in {i,:}(:) ho Ibar della
%                 regione i
%                 totale rimossi, R(t) secondo il SIR
%
%

%[status,result] = fileattrib('dati-andamento-nazionale');
%path_folder = result.Name;                     % percorso alla cartella

path_csv = fullfile(path_folder, '*.csv');
csvfiles=dir(path_csv);                         % struttura contenente campo nome file
numfiles=numel(csvfiles);                       % numero files in path_csv

% inizializzo: creo celle contenenti 21 vettori di 0 del tipo numfiles x 1
regioni = 21;
Ibar = cell(regioni,1);
Rbar = cell(regioni,1);
totCases = cell(regioni,1);
for i = 1:21
    Ibar{i,1} = zeros(numfiles,1);
    Rbar{i,1} = zeros(numfiles,1);
    totCases{i,1}=zeros(numfiles,1);         % numero totale di casi postivi
end

date=string.empty;                      % dichiaro un array vuoto di stringhe

for k=1:numfiles

    path_file = fullfile(path_folder,csvfiles(k).name);  % percorso singolo file dati

	% DATI
    tmp1 = readmatrix(path_file,'Range','K:K');
    tmp2 = sum(readmatrix(path_file,'Range','N:O'),2);   % somma in righe -> ottengo vettore
    tmp3 = readmatrix(path_file,'Range','N:N');     % tot casi positivi

    
    for i = 1:regioni
        Ibar{i,:}(k,:) = tmp1(i);           % totale_positivi = totale_ospedalizzati +
                                                              % isolamento domiciliare
        Rbar{i,:}(k,:) = tmp2(i);           % deceduti + dimessi guariti
        
        totCases{i,:}(k,:) = tmp3(i);       % tot casi positivi

    end
    
    % DATE YYYY/MM/DD in stringhe
    date(k) = readmatrix(path_file,'Range','A2:A2','Delimiter','T','OutputType','string');
    
end

% dall'ultimo file (di cui ho il percorso in path_file) ricavo le regioni
% fid = fopen(path_file,'r');
% fmt  = ['%*s %*s %*d %s %*f %*f' repmat(' %*d',1,13)];
% tmp  = textscan(fid,fmt,'Delimiter',',','HeaderLines',1);
% reg_label = tmp{1};
% fclose(fid);

reg_label = readmatrix(path_file,'Range','D:D','OutputType','string');

% modifico formato date
fmt_date_in = "yyyy-MM-dd"; fmt_date_out = "mmm dd";
DateStringIn = datetime(date,'InputFormat',fmt_date_in);        % converto in datetime
date = datestr(DateStringIn,fmt_date_out);                      % cambio formato data
date = string(date);                                            % casting in string


return
