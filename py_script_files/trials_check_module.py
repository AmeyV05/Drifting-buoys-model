#trial checks for running modules

import readallNC as rdnc



def  main():
 file="../Data_from_models/era5_wind_201403_05.nc"
 [Xa,Ya,Ua,Va,Ta]=rdnc.read_wind(file)
 #[Xt,Yt,Ut,Vt]= [Xt[1344:,1],Yt[1344:,1],Ut[1344:,1],Vt[1344:,1]]
 print(Xa)


if __name__ == '__main__':
	main()