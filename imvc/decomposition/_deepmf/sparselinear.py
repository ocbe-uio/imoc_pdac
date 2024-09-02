import math

try:
    import torch
    import torch.nn
    from torch.nn import Parameter, init
    TORCH_INSTALLED = True
except ImportError:
    class _SparseLinear:
        pass
    TORCH_INSTALLED = False

if TORCH_INSTALLED:
    class _SparseLinear(torch.nn.Module):

        _constants__ = ['bias']

        def __init__(self, in_features, out_features, bias=True):
            super(_SparseLinear, self).__init__()
            self.in_features = in_features
            self.out_features = out_features
            self.weight = Parameter(torch.Tensor(out_features, in_features))
            if bias:
                self.bias = Parameter(torch.Tensor(out_features))
            else:
                self.register_parameter('bias', None)
            self.reset_parameters()

        def reset_parameters(self):
            init.kaiming_uniform_(self.weight, a=math.sqrt(5))
            if self.bias is not None:
                fan_in, _ = init._calculate_fan_in_and_fan_out(self.weight)
                bound = 1 / math.sqrt(fan_in)
                init.uniform_(self.bias, -bound, bound)


        def forward(self, input):

            if input.dim() == 2 and self.bias is not None:
                # fused op is marginally faster
                ret = torch.sparse.addmm(self.bias, input, self.weight.t())
            else:
                output = torch.sparse.mm(input, self.weight.t())
                if self.bias is not None:
                    output += self.bias
                ret = output
            return ret

        def extra_repr(self):
            return 'in_features={}, out_features={}, bias={}'.format(
                self.in_features, self.out_features, self.bias is not None
            )
