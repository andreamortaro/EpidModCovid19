function [date,Ebar,Ibar,Rbar,totCases] = data_read_dpc(path_folder)

%
%   [date,Ibar,Rbar] = data_read_dpc(path_folder)
%
%   data_read_dpc legge i dati della protezione civile italiana
%
%   INPUTS:
%   path_folder : percorso folder dove ci sono i csv
%
%   OUTPUTS:
%   date        : vettore stringhe delle date presenti nei csv
%   Ibar        : totale casi positivi, I(t) secondo il SEIR
%   Ebar        : totale casi positivi, E(t) secondo il SEIR
%   Rbar        : totale rimossi, R(t) secondo il SEIR
%

%[status,result] = fileattrib('dati-andamento-nazionale');
%path_folder = result.Name;         % percorso alla cartella

path_csv = fullfile(path_folder, '*.csv');
csvfiles=dir(path_csv);              % struttura contenente campo nome file
numfiles=numel(dir(path_csv));       % numero files in path_csv

% inizializzo
Ebar=zeros(numfiles,1);
Ibar=zeros(numfiles,1);
Rbar=zeros(numfiles,1);
date=string.empty;                  % dichiaro un array vuoto di stringhe
totCases=zeros(numfiles,1);         % numero totale di casi postivi

for k=1:numfiles

    path_file = fullfile(path_folder,csvfiles(k).name);  % percorso singolo file

    % DATE YYYY/MM/DD in stringhe
    date(k) = readmatrix(path_file,'Range','A2:A2','Delimiter','T','OutputType','string');
    
    % DATI
    Ebar(k) = readmatrix(path_file,'Range','F:F');          % isolamento domiciliare

	Ibar(k) = readmatrix(path_file,'Range','E:E');          % totale_ospedalizzati 
                                                     
    Rbar(k) = sum(readmatrix(path_file,'Range','J:K'));     % deceduti + dimessi guariti
    
    totCases(k) = readmatrix(path_file,'Range','N:N');     % tot casi positivi
        
end

% modifico formato date in vec_date
fmt_date_in = "yyyy-MM-dd"; fmt_date_out = "mmm dd";
DateStringIn = datetime(date,'InputFormat',fmt_date_in);        % converto in datetime
date = datestr(DateStringIn,fmt_date_out);                      % cambio formato data
date = string(date);                                            % casting in string

return
