// Copyright (c) Brandon Li 2025
// This file is part of Brandon's Semiconductor Simulator which is released under GNU GPL v3.0.
// See LICENSE.txt for full license details.

package electrodynamics.util;

import java.util.HashMap;
import java.util.Map;

public class Font7x5 {
    public static final Map<Character, byte[]> FONT_MAP = new HashMap<>();

    static {
        FONT_MAP.put('A', new byte[]{0b01110, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001});
        FONT_MAP.put('B', new byte[]{0b11110, 0b10001, 0b10001, 0b11110, 0b10001, 0b10001, 0b11110});
        FONT_MAP.put('C', new byte[]{0b01110, 0b10001, 0b10000, 0b10000, 0b10000, 0b10001, 0b01110});
        FONT_MAP.put('D', new byte[]{0b11100, 0b10010, 0b10001, 0b10001, 0b10001, 0b10010, 0b11100});
        FONT_MAP.put('E', new byte[]{0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b11111});
        FONT_MAP.put('F', new byte[]{0b11111, 0b10000, 0b10000, 0b11110, 0b10000, 0b10000, 0b10000});
        FONT_MAP.put('G', new byte[]{0b01110, 0b10001, 0b10000, 0b10111, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('H', new byte[]{0b10001, 0b10001, 0b10001, 0b11111, 0b10001, 0b10001, 0b10001});
        FONT_MAP.put('I', new byte[]{0b01110, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110});
        FONT_MAP.put('J', new byte[]{0b00001, 0b00001, 0b00001, 0b00001, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('K', new byte[]{0b10001, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010, 0b10001});
        FONT_MAP.put('L', new byte[]{0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b10000, 0b11111});
        FONT_MAP.put('M', new byte[]{0b10001, 0b11011, 0b10101, 0b10101, 0b10001, 0b10001, 0b10001});
        FONT_MAP.put('N', new byte[]{0b10001, 0b10001, 0b11001, 0b10101, 0b10011, 0b10001, 0b10001});
        FONT_MAP.put('O', new byte[]{0b01110, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('P', new byte[]{0b11110, 0b10001, 0b10001, 0b11110, 0b10000, 0b10000, 0b10000});
        FONT_MAP.put('Q', new byte[]{0b01110, 0b10001, 0b10001, 0b10001, 0b10101, 0b10010, 0b01101});
        FONT_MAP.put('R', new byte[]{0b11110, 0b10001, 0b10001, 0b11110, 0b10100, 0b10010, 0b10001});
        FONT_MAP.put('S', new byte[]{0b01111, 0b10000, 0b10000, 0b01110, 0b00001, 0b00001, 0b11110});
        FONT_MAP.put('T', new byte[]{0b11111, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100});
        FONT_MAP.put('U', new byte[]{0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('V', new byte[]{0b10001, 0b10001, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100});
        FONT_MAP.put('W', new byte[]{0b10001, 0b10001, 0b10001, 0b10101, 0b10101, 0b10101, 0b01010});
        FONT_MAP.put('X', new byte[]{0b10001, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001, 0b10001});
        FONT_MAP.put('Y', new byte[]{0b10001, 0b10001, 0b01010, 0b00100, 0b00100, 0b00100, 0b00100});
        FONT_MAP.put('Z', new byte[]{0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b11111});

        FONT_MAP.put('a', new byte[]{0b00000, 0b00000, 0b01110, 0b00001, 0b01111, 0b10001, 0b01111});
        FONT_MAP.put('b', new byte[]{0b10000, 0b10000, 0b10110, 0b11001, 0b10001, 0b10001, 0b11110});
        FONT_MAP.put('c', new byte[]{0b00000, 0b00000, 0b01110, 0b10001, 0b10000, 0b10001, 0b01110});
        FONT_MAP.put('d', new byte[]{0b00001, 0b00001, 0b01101, 0b10011, 0b10001, 0b10001, 0b01111});
        FONT_MAP.put('e', new byte[]{0b00000, 0b00000, 0b01110, 0b10001, 0b11111, 0b10000, 0b01110});
        FONT_MAP.put('f', new byte[]{0b00110, 0b01001, 0b01000, 0b11110, 0b01000, 0b01000, 0b01000});
        FONT_MAP.put('g', new byte[]{0b00000, 0b00000, 0b01111, 0b10001, 0b01111, 0b00001, 0b11110});
        FONT_MAP.put('h', new byte[]{0b10000, 0b10000, 0b10110, 0b11001, 0b10001, 0b10001, 0b10001});
        FONT_MAP.put('i', new byte[]{0b00100, 0b00000, 0b01100, 0b00100, 0b00100, 0b00100, 0b01110});
        FONT_MAP.put('j', new byte[]{0b00010, 0b00000, 0b00110, 0b00010, 0b00010, 0b10010, 0b01100});
        FONT_MAP.put('k', new byte[]{0b10000, 0b10000, 0b10010, 0b10100, 0b11000, 0b10100, 0b10010});
        FONT_MAP.put('l', new byte[]{0b11000, 0b01000, 0b01000, 0b01000, 0b01000, 0b01001, 0b00110});
        FONT_MAP.put('m', new byte[]{0b00000, 0b00000, 0b11010, 0b10101, 0b10101, 0b10101, 0b10101});
        FONT_MAP.put('n', new byte[]{0b00000, 0b00000, 0b10110, 0b11001, 0b10001, 0b10001, 0b10001});
        FONT_MAP.put('o', new byte[]{0b00000, 0b00000, 0b01110, 0b10001, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('p', new byte[]{0b00000, 0b00000, 0b11110, 0b10001, 0b10001, 0b11110, 0b10000});
        FONT_MAP.put('q', new byte[]{0b00000, 0b00000, 0b01111, 0b10001, 0b10001, 0b01111, 0b00001});
        FONT_MAP.put('r', new byte[]{0b00000, 0b00000, 0b10110, 0b11001, 0b10000, 0b10000, 0b10000});
        FONT_MAP.put('s', new byte[]{0b00000, 0b00000, 0b01111, 0b10000, 0b01110, 0b00001, 0b11110});
        FONT_MAP.put('t', new byte[]{0b01000, 0b01000, 0b11110, 0b01000, 0b01000, 0b01001, 0b00110});
        FONT_MAP.put('u', new byte[]{0b00000, 0b00000, 0b10001, 0b10001, 0b10001, 0b10011, 0b01101});
        FONT_MAP.put('v', new byte[]{0b00000, 0b00000, 0b10001, 0b10001, 0b10001, 0b01010, 0b00100});
        FONT_MAP.put('w', new byte[]{0b00000, 0b00000, 0b10001, 0b10001, 0b10101, 0b10101, 0b01010});
        FONT_MAP.put('x', new byte[]{0b00000, 0b00000, 0b10001, 0b01010, 0b00100, 0b01010, 0b10001});
        FONT_MAP.put('y', new byte[]{0b00000, 0b00000, 0b10001, 0b10001, 0b01111, 0b00001, 0b01110});
        FONT_MAP.put('z', new byte[]{0b00000, 0b00000, 0b11111, 0b00010, 0b00100, 0b01000, 0b11111});

        FONT_MAP.put('0', new byte[]{0b01110, 0b10001, 0b10011, 0b10101, 0b11001, 0b10001, 0b01110});
        FONT_MAP.put('1', new byte[]{0b00100, 0b01100, 0b00100, 0b00100, 0b00100, 0b00100, 0b01110});
        FONT_MAP.put('2', new byte[]{0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b01000, 0b11111});
        FONT_MAP.put('3', new byte[]{0b11111, 0b00010, 0b00100, 0b00010, 0b00001, 0b10001, 0b01110});
        FONT_MAP.put('4', new byte[]{0b00010, 0b00110, 0b01010, 0b10010, 0b11111, 0b00010, 0b00010});
        FONT_MAP.put('5', new byte[]{0b11111, 0b10000, 0b11110, 0b00001, 0b00001, 0b10001, 0b01110});
        FONT_MAP.put('6', new byte[]{0b00110, 0b01000, 0b10000, 0b11110, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('7', new byte[]{0b11111, 0b00001, 0b00010, 0b00100, 0b01000, 0b10000, 0b10000});
        FONT_MAP.put('8', new byte[]{0b01110, 0b10001, 0b10001, 0b01110, 0b10001, 0b10001, 0b01110});
        FONT_MAP.put('9', new byte[]{0b01110, 0b10001, 0b10001, 0b01111, 0b00001, 0b00010, 0b01100});

        FONT_MAP.put(' ', new byte[]{0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000});
        FONT_MAP.put('!', new byte[]{0b00100, 0b00100, 0b00100, 0b00100, 0b00000, 0b00000, 0b00100});
        FONT_MAP.put('"', new byte[]{0b01010, 0b01010, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000});
        FONT_MAP.put('#', new byte[]{0b01010, 0b11111, 0b01010, 0b01010, 0b11111, 0b01010, 0b01010});
        FONT_MAP.put('$', new byte[]{0b00100, 0b01111, 0b10100, 0b01110, 0b00101, 0b11110, 0b00100});
        FONT_MAP.put('%', new byte[]{0b11000, 0b11001, 0b00010, 0b00100, 0b01000, 0b10011, 0b00011});
        FONT_MAP.put('&', new byte[]{0b01000, 0b10100, 0b10100, 0b01000, 0b10101, 0b10010, 0b01101});
        FONT_MAP.put('\'', new byte[]{0b00100, 0b00100, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000});
        FONT_MAP.put('(', new byte[]{0b00010, 0b00100, 0b01000, 0b01000, 0b01000, 0b00100, 0b00010});
        FONT_MAP.put(')', new byte[]{0b01000, 0b00100, 0b00010, 0b00010, 0b00010, 0b00100, 0b01000});
        FONT_MAP.put('*', new byte[]{0b00000, 0b00100, 0b10101, 0b01110, 0b10101, 0b00100, 0b00000});
        FONT_MAP.put('+', new byte[]{0b00000, 0b00100, 0b00100, 0b11111, 0b00100, 0b00100, 0b00000});
        FONT_MAP.put(',', new byte[]{0b00000, 0b00000, 0b00000, 0b00000, 0b00100, 0b00100, 0b01000});
        FONT_MAP.put('-', new byte[]{0b00000, 0b00000, 0b00000, 0b11111, 0b00000, 0b00000, 0b00000});
        FONT_MAP.put('.', new byte[]{0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00100, 0b00100});
        FONT_MAP.put('/', new byte[]{0b00010, 0b00010, 0b00100, 0b00100, 0b00100, 0b01000, 0b01000});
        FONT_MAP.put(':', new byte[]{0b00000, 0b00100, 0b00100, 0b00000, 0b00100, 0b00100, 0b00000});
        FONT_MAP.put(';', new byte[]{0b00000, 0b00100, 0b00100, 0b00000, 0b00100, 0b00100, 0b01000});
        FONT_MAP.put('<', new byte[]{0b00010, 0b00100, 0b01000, 0b10000, 0b01000, 0b00100, 0b00010});
        FONT_MAP.put('=', new byte[]{0b00000, 0b00000, 0b11111, 0b00000, 0b11111, 0b00000, 0b00000});
        FONT_MAP.put('>', new byte[]{0b10000, 0b01000, 0b00100, 0b00010, 0b00100, 0b01000, 0b10000});
        FONT_MAP.put('?', new byte[]{0b01110, 0b10001, 0b00001, 0b00010, 0b00100, 0b00000, 0b00100});
        FONT_MAP.put('@', new byte[]{0b01110, 0b10001, 0b00001, 0b01101, 0b10101, 0b10101, 0b01110});
        FONT_MAP.put('[', new byte[]{0b01110, 0b01000, 0b01000, 0b01000, 0b01000, 0b01000, 0b01110});
        FONT_MAP.put('\\', new byte[]{0b01000, 0b01000, 0b00100, 0b00100, 0b00100, 0b00010, 0b00010});
        FONT_MAP.put(']', new byte[]{0b01110, 0b00010, 0b00010, 0b00010, 0b00010, 0b00010, 0b01110});
        FONT_MAP.put('^', new byte[]{0b00100, 0b01010, 0b10001, 0b00000, 0b00000, 0b00000, 0b00000});
        FONT_MAP.put('_', new byte[]{0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b00000, 0b11111});
        FONT_MAP.put('`', new byte[]{0b01000, 0b00100, 0b00010, 0b00000, 0b00000, 0b00000, 0b00000});
        FONT_MAP.put('{', new byte[]{0b00010, 0b00100, 0b00100, 0b01000, 0b00100, 0b00100, 0b00010});
        FONT_MAP.put('|', new byte[]{0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100, 0b00100});
        FONT_MAP.put('}', new byte[]{0b01000, 0b00100, 0b00100, 0b00010, 0b00100, 0b00100, 0b01000});
        FONT_MAP.put('~', new byte[]{0b00000, 0b00000, 0b00000, 0b01001, 0b10110, 0b00000, 0b00000});
    }

    public static byte[] getCharacter(char c) {
        return FONT_MAP.getOrDefault(c, null);
    }

    public static void printChar(char c) {
        byte[] bitmap = getCharacter(c);
        for (byte row : bitmap) {
            for (int i = 4; i >= 0; i--) {
                System.out.print(((row >> i) & 1) == 1 ? '#' : ' ');
            }
            System.out.println();
        }
    }

    public static int getPixel(char c, int i, int j) {
        byte[] bitmap = getCharacter(c);
        byte row = bitmap[j];
        return (row >> i) & 1;
    }
}