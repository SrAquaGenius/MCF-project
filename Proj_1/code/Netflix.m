ticker = 'nflx.us';

url = ['https://stooq.com/q/d/l/?s=', ticker, '&i=d'];

filename = 'data.csv';
websave(filename, url);

data = readtable(filename);

S = data.Close;
S = S(~isnan(S));