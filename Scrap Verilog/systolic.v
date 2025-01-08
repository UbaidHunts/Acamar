`timescale 1ns / 1ps

module systolic(
    input clk,
    input reset,
    input [31:0] val,
    input [31:0] vec,
    output [31:0]out
    );
  reg [31:0] reg_val;
  reg [31:0] reg_vec;
  wire [31:0] partial_sum_wire;
  reg [31:0] partial_sum;
  
  reg [31:0] product;
  wire [31:0] product_wire;
  
  fp_multiplier fpm(
     .a(reg_val),
     .b(reg_vec),
     .result(product_wire)
     );
     
   fp_adder fpa(
     .a(product),
     .b(partial_sum),
     .result(partial_sum_wire)
     );

   always @(posedge clk or posedge reset) begin
        if (reset) begin
            reg_val <= 32'b0;
            reg_vec <= 32'b0;
            product <= 32'b0;
            partial_sum <= 32'b0;
        end else begin
            reg_val <= val;
            reg_vec <= vec;
            product<=product_wire;
            partial_sum<=partial_sum_wire;
        end
    end 
    assign out = reg_vec;      
endmodule

module fp_adder (
    input wire [31:0] a, // First floating-point operand
    input wire [31:0] b, // Second floating-point operand
    output reg [31:0] result // Result of the addition
);

    reg sign_a, sign_b, sign_res;
    reg [7:0] exponent_a, exponent_b, exponent_res;
    reg [23:0] mantissa_a, mantissa_b, mantissa_res;
    reg [24:0] mantissa_sum;
    reg [7:0] exponent_diff;
    reg [23:0] mantissa_tmp;
    reg carry;

    always @(*) begin
        // Extract sign, exponent, and mantissa
        sign_a = a[31];
        sign_b = b[31];
        exponent_a = a[30:23];
        exponent_b = b[30:23];
        mantissa_a = {1'b1, a[22:0]}; // Add implicit leading 1
        mantissa_b = {1'b1, b[22:0]}; // Add implicit leading 1

        // Align exponents
        if (exponent_a > exponent_b) begin
            exponent_diff = exponent_a - exponent_b;
            mantissa_b = mantissa_b >> exponent_diff;
            exponent_res = exponent_a;
        end else begin
            exponent_diff = exponent_b - exponent_a;
            mantissa_a = mantissa_a >> exponent_diff;
            exponent_res = exponent_b;
        end

        // Add or subtract mantissas
        if (sign_a == sign_b) begin
            mantissa_sum = mantissa_a + mantissa_b;
            sign_res = sign_a;
            carry = mantissa_sum[24]; // Check for carry out
            if (carry) begin
                mantissa_res = mantissa_sum[24:1]; // Normalize the mantissa
                exponent_res = exponent_res + 1;
            end else begin
                mantissa_res = mantissa_sum[23:0];
            end
        end else begin
            if (mantissa_a > mantissa_b) begin
                mantissa_tmp = mantissa_a - mantissa_b;
                sign_res = sign_a;
            end else begin
                mantissa_tmp = mantissa_b - mantissa_a;
                sign_res = sign_b;
            end
            // Normalize the result if necessary
            if (mantissa_tmp[23]) begin
                mantissa_res = mantissa_tmp;
            end else begin
                mantissa_res = mantissa_tmp << 1;
                exponent_res = exponent_res - 1;
            end
        end

        // Pack the result
        result = {sign_res, exponent_res, mantissa_res[22:0]};
    end
endmodule

module fp_multiplier (
    input wire [31:0] a,
    input wire [31:0] b,
    output reg [31:0] result
);

    reg sign_a;
    reg sign_b;
    reg [7:0] exponent_a;
    reg [7:0] exponent_b;
    reg [23:0] mantissa_a;
    reg [23:0] mantissa_b;
    reg sign_result;
    reg [8:0] exponent_result;
    reg [47:0] mantissa_product;

    always @(*) begin
        // Extract sign, exponent, and mantissa
        sign_a = a[31];
        sign_b = b[31];
        exponent_a = a[30:23];
        exponent_b = b[30:23];
        mantissa_a = {1'b1, a[22:0]}; // Add implicit leading 1
        mantissa_b = {1'b1, b[22:0]}; // Add implicit leading 1

        // Compute the result sign
        sign_result = sign_a ^ sign_b;

        // Compute the result exponent
        exponent_result = exponent_a + exponent_b - 127;

        // Multiply the mantissas
        mantissa_product = mantissa_a * mantissa_b;

        // Normalize the result
        if (mantissa_product[47]) begin
            exponent_result = exponent_result + 1;
            mantissa_product = mantissa_product >> 1;
        end

        // Pack the result
        result = {sign_result, exponent_result[7:0], mantissa_product[46:24]};
    end
endmodule