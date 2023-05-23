mises_ = (pow(sigma_.component(0)-sigma_.component(3),2)+pow(sigma_.component(0)-sigma_.component(5),2)+pow(sigma_.component(3)-sigma_.component(5),2))/2+3*(pow(sigma_.component(1),2)+pow(sigma_.component(2),2)+pow(sigma_.component(4),2));

infovstime.open ("infoSolid.csv", std::ios::out|std::ios::trunc);
infovstime << "%Time,minDx,minDy,minDz,maxDx,maxDy,maxDz,minsigma_xx,minsigma_xy,minsigmayy,maxsigma_xx,maxsigma_xy,maxsigma_yy,maxmises" << nl;
infovstime << runTime().timeName() << ","
           << min(D_.component(0)).value() << ","
           << min(D_.component(1)).value() << ","
           << min(D_.component(2)).value() << ","
           << max(D_.component(0)).value() << ","
           << max(D_.component(1)).value() << ","
           << max(D_.component(2)).value() << ","
           << min(sigma_.component(0)).value() << ","
           << min(sigma_.component(1)).value() << ","
           << min(sigma_.component(3)).value() << ","
           << max(sigma_.component(0)).value() << ","
           << max(sigma_.component(1)).value() << ","
           << max(sigma_.component(3)).value() << ","
           << pow(max(mises_).value(),0.5) << nl;
           
infovstime.flush();
