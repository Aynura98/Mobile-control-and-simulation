classdef EKF < handle
    
    properties (Access=public)
        Xh_kk1;
        Xh_kk;
        P_kk1;
        P_kk;
        
        V;
        W;
        
        
        Xh_kk1_vec;
        Xh_kk_vec;
        P_kk1_vec;
        P_kk_vec;
        Kk_vec;
        innovation_vec;
        
        

        f_fun;
        h_fun;
        Ak_fun;
        Ck_fun;
        
    end
    
    
    methods (Access=public)
        function[obj]=EKF(f_fun,h_fun,Ak_fun,Ck_fun,Xhat0,P0,V,W)
           
            obj.f_fun=f_fun;
            obj.h_fun=h_fun;
            obj.Ak_fun=Ak_fun;
            obj.Ck_fun=Ck_fun;
            
            %% filter matrices
            obj.V=V;
            obj.W=W;
            obj.P_kk1=P0;
            
            %% filter history
            obj.Xh_kk1=Xhat0;
            obj.Xh_kk1_vec=Xhat0;
            obj.P_kk1_vec=P0;
            obj.P_kk_vec=[];
            

        end
 
        function[X_hat_kk,P_kk,innovation,Kk]=estimate(obj,yk,uk)
            Xhkk1=obj.Xh_kk1; 
            Pkk1=obj.P_kk1;
            VV=obj.V;
            
            Ck=obj.Ck_fun(Xhkk1,uk);%template method

            
            %% kalman gain
            Kk=Pkk1*Ck'/((Ck*Pkk1*Ck'+VV));

            %% prefdicted output
            yhk=obj.h_fun(Xhkk1,uk);
            
            %% innovation
            ee=yk-yhk;
            
            %% from prediction to estimation
            Xhkk=Xhkk1+Kk*ee;
            Pkk=Pkk1-Kk*Ck*Pkk1;
            
            %% output
            innovation=ee;
            X_hat_kk=Xhkk;
            P_kk=Pkk;
            
            %% storing data
            obj.Xh_kk=Xhkk;
            obj.Xh_kk_vec=[obj.Xh_kk_vec Xhkk];
            obj.P_kk=Pkk;
            if isempty(obj.P_kk_vec)
                N=0;
            else
                N=size(obj.P_kk_vec,3);
            end
            obj.P_kk_vec(:,:,N+1)=Pkk;
            obj.Kk_vec(:,:,N+1)=Kk;
            obj.innovation_vec=[obj.innovation_vec ee];
        end
        
        function[X_hat_k1k,P_k1k]=predict(obj,uk)
            
            Xhkk=obj.Xh_kk;
            Pkk=obj.P_kk;
            WW=obj.W;
            
            Ak=obj.Ak_fun(Xhkk,uk);
            
            
            %% from estimation to prediction
            Xhk1k=obj.f_fun(Xhkk,uk);
            Pk1k=Ak*Pkk*Ak'+WW;
            
            %% output
            X_hat_k1k=Xhk1k;
            P_k1k=Pk1k;
            
            %% stroring data
            obj.Xh_kk1=Xhk1k;
            obj.Xh_kk1_vec=[obj.Xh_kk1_vec Xhk1k];
            obj.P_kk1=Pk1k;
            if isempty(obj.P_kk1_vec)
                N=0;
            else
                N=size(obj.P_kk1_vec,3);
            end
            obj.P_kk1_vec(:,:,N+1)=Pk1k;
        end
        
        
        function[X_hat_kk,P_kk]=update(obj,yk,uk)
            %% estimation
            [X_hat_kk,P_kk,~,~]=obj.estimate(yk,uk);
            
            %% prediction
            obj.predict(uk);
        end
        
        function[X_hat_kk_vec,P_kk_vec]=get_kk_data(obj)
            X_hat_kk_vec=obj.Xh_kk_vec;
            P_kk_vec=obj.P_kk_vec;
        end
        function[X_hat_k1k_vec,P_k1k_vec]=get_k1k_data(obj)
            X_hat_k1k_vec=obj.Xh_kk1_vec;
            P_k1k_vec=obj.P_kk1_vec;
        end
        function[innovation_vec]=get_innovation(obj)
            innovation_vec=obj.innovation_vec;
        end
        function[Kk_vec]=get_Kalman_gains(obj)
            Kk_vec=obj.Kk_vec;
        end
        
    end
    methods (Access=private)
    
    end
    
end

