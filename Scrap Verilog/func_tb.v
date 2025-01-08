`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: 
// Engineer: 
// 
// Create Date: 06/14/2024 04:08:59 PM
// Design Name: 
// Module Name: func_tb
// Project Name: 
// Target Devices: 
// Tool Versions: 
// Description: 
// 
// Dependencies: 
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////

/*module func_tb();

    // Parameters
    parameter CLK_PERIOD = 10;  // Clock period in time units

    // Inputs
    reg clk;
    reg reset;
    reg [31:0] t_val;
    reg [31:0] t_vec;

    // Outputs
    wire [31:0] t_out;

    // Instantiate the PE module
    systolic uut (
        .clk(clk),
        .reset(reset),
        .val(t_val),
        .vec(t_vec),
        .out(t_out)
    );

    // Clock generation
    always #(CLK_PERIOD / 2) clk = ~clk;

    // Test procedure
    initial begin
        // Initialize signals
        clk = 0;
        reset = 1;
        t_val = 32'b0;
        t_vec = 32'b0;

        // Apply reset
        #(CLK_PERIOD * 2);
        reset = 0;

        // Test vector 1
        t_val = 32'd1;  // 1.0 in IEEE 754 single precision
        t_vec = 32'd5;  // 2.0 in IEEE 754 single precision
        #(CLK_PERIOD * 2);
        

        // Test vector 2
        t_val = 32'd0; // 1.0 in IEEE 754 single precision
        t_vec = 32'd3;  // 1.0 in IEEE 754 single precision
        #(CLK_PERIOD * 2);


        // Test vector 3
        t_val = 32'd3; // 0.5 in IEEE 754 single precision
        t_vec = 32'd4;  // 3.0 in IEEE 754 single precision
        #(CLK_PERIOD * 2);
product

        // End of simulation
        $stop;
    end

endmodule*/

module testbench;

    // Parameters
    parameter CLK_PERIOD = 10; // Clock period in ns

    // Signals
    reg clk;
    reg reset;
    reg [31:0] val1;
    reg [31:0] val2;
    reg [11:0] rowIdx1;
    reg [11:0] rowIdx2;
    reg tag1;
    reg tag2;
    reg [31:0] vec;
    wire overlap;

    // Instantiate the module under test
    conf_sys dut (
        .clk(clk),
        .reset(reset),
        .val1(val1),
        .val2(val2),
        .rowIdx1(rowIdx1),
        .rowIdx2(rowIdx2),
        .tag1(tag1),
        .tag2(tag2),
        .vec(vec),
        .overlap(overlap)
    );

    // Clock generation
    always #(CLK_PERIOD / 2) clk = ~clk;


    // Initial stimulus
    initial begin
        // Reset sequence
        reset = 1;
        val1 = 32'b0;
        val2 = 32'b0;
        rowIdx1 = 12'b0;
        rowIdx2 = 12'b0;
        tag1 = 1'b0;
        tag2 = 1'b1;
        vec = 32'b0;
        reset = 0;
        #(CLK_PERIOD * 2);
        
        // Test case 1: val1 and val2 non-zero
        val1 = 32'hABCDE;  // Example non-zero values
        val2 = 32'h54321;
        rowIdx1 = 12'h123;
        rowIdx2 = 12'h456;
        tag1 = 1'b1;
        tag2 = 1'b0;
        vec = 32'h1234;
        #(CLK_PERIOD * 2);

        // Test case 2: val1 zero, val2 non-zero
        val1 = 32'b0;
        val2 = 32'h98765;
        rowIdx1 = 12'b0;
        rowIdx2 = 12'h789;
        tag1 = 1'b0;
        tag2 = 1'b1;
        vec = 32'h5678;
        #(CLK_PERIOD * 2);
        // Add more test cases as needed

        // End simulation
        $finish;
    end

endmodule
