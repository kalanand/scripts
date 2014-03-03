        subroutine getinfo(h1, h2)
        real prob
        vector getprob(1)
        call hdiff(h1, h2, prob)
          getprob(1)=prob
        end
