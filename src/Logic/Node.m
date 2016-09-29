classdef (Abstract,HandleCompatible) Node < handle & matlab.mixin.Heterogeneous
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        children;
        parent;        
    end
    
    methods(Abstract)        
        evaluate(self,assignment);        
        res = toString(self,PipeAnd);
        dnfNode = convertToDNF(self);    
        reduce(self);
%         contains(self,literal)        
    end
    methods
        function obj = Node()
            obj.children = [];
        end
        
        function addChild(self,childNode)
            if isempty(self.children)
                self.children = childNode;
            else
                self.children(end+1) = childNode;
            end            
            childNode.parent = self;
        end                
        function res = contains(self,literal)
            for c=1:numel(self.children)
                child = self.children(c);
                res = child.contains(literal);
                if res
                    break;
                end
            end
        end
        
        function literals = getLiterals(self) 
            literals = [];
            for c=1:numel(self.children)
                child = self.children(c);
                literals = [literals child.getLiterals()];
            end
        end
    end
    
end

