function [ positions ] = subplot_pos(plotwidth,plotheight,leftmargin,rightmargin,bottommargin,topmargin,ncol,nrow,spacex,spacey)
 
    subxsize=(plotwidth-leftmargin-rightmargin-spacex*(ncol-1.0))/ncol;
    subysize=(plotheight-topmargin-bottommargin-spacey*(nrow-1.0))/nrow;
 
    
    for j=1:nrow
        for i=1:ncol
            xfirst=leftmargin+(i-1.0)*(subxsize+spacex);
            yfirst=bottommargin+(j-1.0)*(subysize+spacey);
            
            positions{j,i}=[xfirst/plotwidth yfirst/plotheight subxsize/plotwidth subysize/plotheight];
            
        end
    end
end

